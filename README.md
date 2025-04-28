# Metabarcoding

Unified eDNA Metabarcoding Pipeline (Nanopore) – Consolidated Bash Script

This section outlines a streamlined single Bash script that integrates the user’s eDNA metabarcoding workflow with improvements inspired by the PIMENTA pipeline and recent protocols (e.g., Van der Vorst et al. 2024, Erkenswick et al. 2024). The pipeline is organized into clearly labeled phases and implements best practices (skipping completed steps, robust Bash syntax for WSL compatibility, etc.) for ease of use and reproducibility.
Phase 1: Read Preprocessing and Quality Filtering
In Phase 1, raw Nanopore reads (already basecalled and demultiplexed by Guppy) are preprocessed. This includes poly-A/T tail removal and stringent quality/length filtering, while deferring any primer trimming to a later stage:
Input: Demultiplexed FASTQ files per sample (Guppy has trimmed adapters during basecalling). Reads from all targeted loci remain combined at this stage.


Poly-A/T Trimming: Use Cutadapt to trim homopolymer tails (common in Nanopore reads) – e.g., remove runs of AAAA... or TTTT... at read ends ​(github.com). This step parallels PIMENTA’s initial adapter and poly-A trimming (PIMENTA uses Guppy + PRINSEQ for similar effect (​github.com).


Quality Filtering: Apply minimum quality thresholds (e.g., mean Q-score cutoff) and length range filters. For example, require reads within an expected amplicon length range and discard low-quality reads (Cutadapt’s -q and -m/-M options can enforce this). In PIMENTA, PRINSEQ performed an analogous trimming of low-quality bases and enforced length bounds (​github.com).


Output: Cleaned FASTQ files per sample (postfix e.g. _trimmed.fastq). These files contain full-length amplicon reads with adapters removed and homopolymers trimmed, but still include primer sequences at their ends (primer removal is postponed to Phase 3).


Note: Deferring primer trimming until after consensus generation helps preserve maximum overlap during clustering/alignment. This ensures that even reads of different marker loci remain separated by inherent sequence dissimilarity rather than by manual partitioning, and maximizes the data retained for consensus building ​(github.com).
Each sample’s preprocessing can be wrapped in a check – for example, the script can skip this phase for a sample if its _trimmed.fastq already exists (use if [ -f "$sample_trimmed.fastq" ] guard). This way, reruns don’t redo Phase 1 unnecessarily.
Phase 2: Per-Sample Clustering & Consensus (NGSpeciesID)
In Phase 2, each sample’s filtered reads are clustered by sequence similarity and collapsed into high-accuracy consensus sequences per barcode (marker) using NGSpeciesID. This tool implements reference-free clustering (e.g., via isONclust or a similar algorithm) followed by multiple sequence alignment and consensus calling, analogous to the first clustering step of PIMENTA (​github.com):
Clustering per Sample: NGSpeciesID groups reads into clusters based on sequence similarity (e.g., reads ≥90% similar might form a cluster). This is similar to other long-read metabarcoding approaches which cluster reads at ~90% identity to delineate putative Molecular OTUs​ (bmcgenomics.biomedcentral.com). Each cluster ideally represents a unique DNA barcode (species or taxon) in that sample.


Consensus Calling: For each cluster, NGSpeciesID builds a multiple alignment (using an internal aligner like SPOA) and calls a majority-rule consensus sequence. This consensus corrects random ONT errors, yielding sequences with accuracy often >99% (​bmcgenomics.biomedcentral.com). Recent studies show that consensus barcodes from Nanopore reads can be as accurate and indel-free as Illumina sequences when such clustering and consensus algorithms are applied (​bmcgenomics.biomedcentral.com​bmcgenomics.biomedcentral.com).


Multiple Consensus per Sample: NGSpeciesID may output multiple consensus sequences per sample, one per distinct cluster/taxon detected (​bmcbiol.biomedcentral.com​bmcbiol.biomedcentral.com). This is crucial for eDNA/metabarcoding, as a single sample (environmental sample) can contain DNA from many species. (For each input “barcode”/sample, NGSpeciesID will produce 0, 1, or many consensus sequences depending on diversity.)


Output: For each sample and each primer-defined locus, a FASTA file of consensus sequences (e.g., sampleX_consensus.fasta). Each FASTA header can carry the sample ID and cluster ID. Example: >SampleX_cluster1 readCount=50 followed by the consensus sequence.
 
NGSpeciesID automates clustering and consensus in one step (​bmcbiol.biomedcentral.com), so it replaces the manual CD-HIT → MAFFT → consensus steps that some pipelines use for initial clustering. If an expected per-sample consensus FASTA already exists, the script will skip rerunning NGSpeciesID for that sample (preventing redundant computation).
Phase 3: Cross-Sample Consensus Deduplication (Final Clustering & Primer Trimming)
Phase 3 takes all consensus sequences from Phase 2 (pooled across samples, for all barcode loci) and performs a final clustering and deduplication, followed by primer trimming. This step collapses identical or nearly identical sequences across samples and ensures each unique sequence (presumably representing a unique taxon DNA barcode) is represented once. We use a two-step approach inspired by the PIMENTA “reclustering” phase​ (github.com):
All-vs-All Clustering (CD-HIT): We combine all per-sample consensus FASTAs into one file and run CD-HIT-EST at a high identity threshold (e.g., 98–100%). This clusters redundant consensus sequences across samples (github.com). For instance, if the same fish DNA sequence was present in multiple samples, CD-HIT will group them. We choose a stringent cutoff (near 100%) so that only almost-identical sequences cluster, preserving true biological variants as separate clusters​ (bmcgenomics.biomedcentral.com). (PIMENTA performs an initial 90% pre-cluster followed by a second high-identity clustering (​github.com), but here a single high-identity pass is sufficient for deduplication.)


Consensus per Cluster: For each cluster produced by CD-HIT, we align the member sequences with MAFFT and build a final consensus sequence using the user’s consensus.pl script (similar to PIMENTA’s consensus generation step)​ (github.com). This yields one representative sequence per cluster (ensuring any minor discrepancies between samples are resolved). If CD-HIT was set to 100% identity, the consensus is trivial (identical sequences); if 98–99%, the consensus will correct any single-base differences among nearly identical sequences.


Primer Removal (Cutadapt): Now, trim primer sequences from these final consensus sequences. We utilize Cutadapt only at this stage to remove the known primer binding regions from the consensus ends. Multiple primer pairs can be handled by grouping sequences by primer set:


Define a Bash array PRIMER_PAIRS listing each locus and its forward/reverse primer sequences (e.g., PRIMER_PAIRS=( "Fish16S AGCTCCTGAAA... TGGTGCAGTT..." "COI GGTCAACAAATC... CATCCAACAT..." )).


For each primer pair, run Cutadapt on the combined consensus FASTA, specifying the forward primer as a 5′ adapter and the reverse primer as a 3′ adapter (using -g/--front and -a/--adapter with appropriate anchoring). Use options to allow partial matches and strip them from sequences. Important: ensure Cutadapt outputs trimmed sequences for that primer set (e.g., using -o marker_trimmed.fasta and perhaps --untrimmed-output for those that didn’t have that primer).


Repeat for each primer pair, then concatenate or otherwise organize the trimmed outputs. Each consensus sequence should end up in the file corresponding to its correct primer set (since non-matching sequences go to the “untrimmed” file in each round, they’ll be captured when their correct primer pair is processed). This effectively sorts final consensus sequences by marker and removes primer regions.


Output: Final primer-trimmed consensus FASTA files (or one combined FASTA) containing unique high-quality barcode sequences ready for taxonomic assignment. Primer sequences (which can confound BLAST by matching many reference entries) are now removed, per the requirement to trim primers only post-consensus. By iterating over PRIMER_PAIRS, the script cleanly handles multi-locus datasets.


Throughout Phase 3, the script uses clear subdirectories (e.g., Clustering/All_consensus.fasta, Clustering/CDHIT_clusters.clstr, Consensus/final_consensus.fasta, Consensus/trimmed_by_marker/) to organize intermediate files. Skipping logic can be applied at each sub-step; for example, if final_consensus.fasta already exists from a previous run, the CD-HIT and MAFFT steps can be skipped.
Phase 4: Taxonomic Assignment (BLAST + Post-processing)****
In Phase 4, each unique consensus sequence is subjected to taxonomic identification using a BLAST search against NCBI’s nucleotide database, followed by filtering of results to ensure reliable assignments. The BLAST and taxonomic classification procedure includes:
BLAST (megablast) Search: Use blastn in megablast mode (optimized for highly similar sequences) to query each consensus sequence against the NCBI nt database. We utilize parameters to limit results and speed up searches, e.g., -max_target_seqs 10 (retrieve top 10 hits) and an e-value cutoff (e.g., 1e-5 or stricter). The pipeline focuses on near-exact matches by design, so megablast is appropriate (​github.com).


Identity and Coverage Filters: Parse the BLAST output to retain only hits with ≥90% sequence identity AND ≥90% query coverage​( github.com). This filter (PIMENTA’s default is 90/90 as well) ensures we only consider hits that match almost the entire query at high similarity – a common threshold in eDNA for confident species-level matches. Any hit below either threshold is ignored in downstream analysis.


Exclude Unwanted Taxa: Implement a filter to remove hits that correspond to environmental or unplaced sequences, using user-provided exclusion lists (e.g., a list of Taxonomy IDs or keywords for contaminants). For example, any hit with a taxonomy containing terms like “uncultured”, “environmental sample”, or matching TaxIDs in an exclusion file is discarded. (The PIMENTA pipeline provides Excluded.NCBI.identifications.tsv and Excluded.NCBI.taxids.tsv to facilitate this filtering (​github.com).) This prevents assignments to ambiguous environmental entries in GenBank.


Consensus Taxonomy & Rank Adjustment: For each query, determine the consensus taxonomy from the remaining top hits:


If the top hits (after filtering) all correspond to the same species, assign that species with confidence.


If multiple top hits are present for different species within the same genus, downgrade to genus-level identification. This scenario indicates ambiguity at the species level (the sequence could be one of several closely related species) – the script will report the genus instead of picking a false species match.


If top hits span multiple genera (or higher taxa), further reduce the taxonomic resolution (e.g., report family or order, or mark as “ambiguous”). Essentially, the script finds the lowest common taxonomic rank shared by all top hits above the identity/coverage threshold. For instance, if one top hit is Species A (Genus X) and another is Species B (Genus Y), but both are in the same Family, the assignment can be given at Family level; if they differ at family, assign at Order, etc. This dynamic rank adjustment ensures we report conservative identifications when BLAST cannot unequivocally pinpoint a single species (​github.com).


Output: A tabular report (and/or annotated FASTA) listing each consensus sequence and its assigned taxonomy (with taxonomic rank noted), plus BLAST metrics (identity, coverage, accession of top hit). The script can also produce a human-readable summary, and even Krona charts or similar (as PIMENTA does with Krona for visualization (​github.com)).


This BLAST-based taxonomy step is designed to maximize accuracy. By requiring high percent identity and coverage and by excluding problematic reference sequences, we follow the strict practices recommended in recent eDNA metabarcoding studies for minimizing false positives (github.com). The logic for taxonomic resolution mirrors approaches in the literature that emphasize not over-calling species when sequences are conserved across taxa (common in short barcodes).
Finally, the script can incorporate a skip check for BLAST: if a previous run’s results file (e.g., blast_results.tsv or a final taxonomy report) exists for a given consensus FASTA, the BLAST search will be skipped or results merged, to avoid re-querying the database for identical sequences.
