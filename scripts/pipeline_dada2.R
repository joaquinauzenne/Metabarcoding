#!/usr/bin/env Rscript
# pipeline_dada2.R — DADA2 ASV inference from Nanopore reads
#
# EXPERIMENTAL: DADA2's error model was designed for Illumina data.
# This script tests whether DADA2 can learn a useful error model from
# Nanopore reads that have been quality-filtered but NOT consensus-polished.
#
# After running, inspect error_rates.pdf:
#   - Clear decreasing trend with quality score → model is usable
#   - Flat or chaotic rates → DADA2 hasn't learned the error structure;
#     treat its ASVs as unreliable and defer to NGSpeciesID-ASV results
#
# Usage: Rscript pipeline_dada2.R <trimmed_dir> <output_dir> <threads>
#
# Input:  Quality-filtered FASTQ files from Phase 1 (one per barcode)
#         (i.e., *_trimmed.fastq.gz files from the Cutadapt step)
# Output: dada2_asvs.fasta         — ASV sequences for BLAST
#         dada2_asv_table.csv      — samples × ASVs abundance matrix
#         dada2_summary.csv        — per-sample read counts
#         quality_profiles.pdf     — Q-score profiles (first 4 samples)
#         error_rates.pdf          — learned error rates (INSPECT THIS FIRST)

suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript pipeline_dada2.R <trimmed_dir> <output_dir> [threads]",
       call. = FALSE)
}
trimmed_dir <- args[1]
output_dir  <- args[2]
threads     <- if (length(args) >= 3) as.integer(args[3]) else 1L

# ── Primer set config ──────────────────────────────────────────────────────────
# Set PRIMER_SET, MIN_LEN, MAX_LEN, and MAX_EE to match your library.
# Run the script, then inspect quality_profiles.pdf to validate MAX_EE.
#
# Supported primer sets and expected amplicon sizes:
#   12S — Ushio et al. 2017:        ~160 bp → MIN_LEN = 130, MAX_LEN = 200
#   12S — Taberlet et al. 2018:      ~60 bp → MIN_LEN = 45,  MAX_LEN = 90
#   12S — Riaz et al. 2011:         ~107 bp → MIN_LEN = 85,  MAX_LEN = 140
#   16S — Boessenkool et al. 2012:   ~90 bp → MIN_LEN = 70,  MAX_LEN = 120
#   16S — Lynggaard et al. 2024:    ~300 bp → MIN_LEN = 250, MAX_LEN = 370
#   16S — Franklin et al. 2019:     ~330 bp → MIN_LEN = 280, MAX_LEN = 400
#
# *** SET THESE BEFORE RUNNING — do not leave at default from a previous run ***
PRIMER_SET <- "16S_Boessenkool"  # change to match your actual library
MIN_LEN    <- 70L                # lower bound (bp) — ~80% of expected amplicon
MAX_LEN    <- 120L               # upper bound (bp) — ~130% of expected amplicon

# maxEE: calibrate against NanoPlot quality distributions from Phase 2.
# At Q12 (typical Nanopore), expected errors accumulate ~1 per 16 bp.
# Verify against quality_profiles.pdf before committing to this value.
MAX_EE <- 5

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ── Step 1: Identify input files ───────────────────────────────────────────────
# Nanopore = single-end reads (no paired-end merging)
fastq_files <- sort(list.files(trimmed_dir,
                                pattern = "\\.fastq$|\\.fastq\\.gz$",
                                full.names = TRUE,
                                recursive = TRUE))

if (length(fastq_files) == 0) {
  stop("No FASTQ files found in: ", trimmed_dir)
}

sample_names <- sub("_trimmed\\.fastq.*", "", basename(fastq_files))
names(fastq_files) <- sample_names
message(sprintf("Found %d samples", length(fastq_files)))

# ── Step 2: Quality profiles (inspect before committing to MAX_EE) ────────────
pdf(file.path(output_dir, "quality_profiles.pdf"))
n_plot <- min(4L, length(fastq_files))
for (fn in fastq_files[seq_len(n_plot)]) {
  print(plotQualityProfile(fn))
}
dev.off()
message("Quality profiles written to: ",
        file.path(output_dir, "quality_profiles.pdf"))

# ── Step 3: Filter and trim ────────────────────────────────────────────────────
# maxN = 0: DADA2's error model cannot handle ambiguous bases
# truncQ = 2: trim at the first Q≤2 base (noise-only tail removal)
# rm.phix = FALSE: PhiX is an Illumina spike-in; irrelevant for Nanopore
filt_dir   <- file.path(output_dir, "filtered")
filt_files <- file.path(filt_dir, paste0(sample_names, "_filt.fastq.gz"))
names(filt_files) <- sample_names

message(sprintf("Filtering: PRIMER_SET=%s | minLen=%d | maxLen=%d | maxEE=%g",
                PRIMER_SET, MIN_LEN, MAX_LEN, MAX_EE))

filt_out <- filterAndTrim(
  fastq_files, filt_files,
  minLen      = MIN_LEN,
  maxLen      = MAX_LEN,
  maxN        = 0,
  maxEE       = MAX_EE,
  truncQ      = 2,
  rm.phix     = FALSE,
  multithread = threads,
  verbose     = TRUE
)

exists <- file.exists(filt_files)
if (sum(exists) == 0) {
  stop("All samples failed filtering. Check MIN_LEN/MAX_LEN/MAX_EE settings.")
}
filt_files   <- filt_files[exists]
sample_names <- sample_names[exists]
message(sprintf("%d / %d samples passed filtering", sum(exists), length(exists)))

# ── Step 4: Learn error rates ──────────────────────────────────────────────────
# KEY STEP for Nanopore compatibility. The model may not fit as cleanly as for
# Illumina but should still show a decreasing trend with quality score.
message("Learning error rates (this may take a while)...")
err <- learnErrors(filt_files,
                   multithread  = threads,
                   randomize    = TRUE,   # sample reads randomly
                   MAX_CONSIST  = 20,     # allow more iterations to converge
                   verbose      = TRUE)

# Save error plot — compare to Illumina DADA2 tutorial plots.
# Nanopore plots will look different (higher overall error rate, especially
# for indels) but should show a decreasing trend with quality score.
pdf(file.path(output_dir, "error_rates.pdf"))
print(plotErrors(err, nominalQ = TRUE))
dev.off()
message("Error rates written to: ", file.path(output_dir, "error_rates.pdf"))
message("*** INSPECT error_rates.pdf BEFORE TRUSTING RESULTS ***")

# ── Step 5: Dereplicate ────────────────────────────────────────────────────────
derep <- derepFastq(filt_files, verbose = TRUE)
if (length(filt_files) == 1L) {
  derep <- list(derep)
  names(derep) <- sample_names
}

# ── Step 6: Sample inference (core DADA2 algorithm) ───────────────────────────
# pool = TRUE: cross-sample pooling helps detect low-abundance ASVs,
# important for rare species in eDNA
message("Running DADA2 sample inference (pool=TRUE)...")
dada_out <- dada(derep, err = err,
                 multithread = threads,
                 pool        = TRUE,
                 verbose     = TRUE)

# ── Step 7: Construct ASV table ────────────────────────────────────────────────
seqtab <- makeSequenceTable(dada_out)
message(sprintf("ASV table: %d samples × %d ASVs",
                nrow(seqtab), ncol(seqtab)))

# ── Step 8: Remove chimeras ────────────────────────────────────────────────────
seqtab_nochim <- removeBimeraDenovo(seqtab,
                                     method      = "consensus",
                                     multithread = threads,
                                     verbose     = TRUE)
pct_retained <- sum(seqtab_nochim) / sum(seqtab) * 100
message(sprintf("%.1f%% of reads retained after chimera removal", pct_retained))
message(sprintf("Final: %d ASVs", ncol(seqtab_nochim)))

# ── Step 9: Write outputs ──────────────────────────────────────────────────────
asv_seqs <- colnames(seqtab_nochim)
asv_ids  <- paste0("dada2_asv_", seq_along(asv_seqs))

# ASV sequences as FASTA (for BLAST in pipeline.sh Phase 4)
fasta_out <- file.path(output_dir, "dada2_asvs.fasta")
writeLines(paste0(">", asv_ids, "\n", asv_seqs), fasta_out)
message("ASV FASTA written to: ", fasta_out)

# ASV abundance table (samples × ASVs)
asv_table           <- as.data.frame(seqtab_nochim)
colnames(asv_table) <- asv_ids
asv_table$sample    <- rownames(asv_table)
asv_table_out <- file.path(output_dir, "dada2_asv_table.csv")
write.csv(asv_table, asv_table_out, row.names = FALSE)
message("ASV table written to: ", asv_table_out)

# Summary statistics
summary_df <- data.frame(
  sample         = sample_names,
  input_reads    = filt_out[exists, "reads.in"],
  filtered_reads = filt_out[exists, "reads.out"],
  denoised_reads = sapply(dada_out, function(x) sum(x$denoised)),
  n_asvs         = rowSums(seqtab_nochim > 0)
)
summary_out <- file.path(output_dir, "dada2_summary.csv")
write.csv(summary_df, summary_out, row.names = FALSE)
message("Summary written to: ", summary_out)

message(sprintf(
  "\nDADA2 complete. %d ASVs across %d samples. Inspect error_rates.pdf before trusting results.",
  length(asv_seqs), nrow(seqtab_nochim)
))
message(sprintf(
  "Next step: run BLAST on %s, then Rscript scripts/lca_assignment.R", fasta_out
))
