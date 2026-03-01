#!/bin/bash
#SBATCH --job-name=mammalia_nt_db       # Job name
#SBATCH --output=mammalia_nt.%j.out     # STDOUT
#SBATCH --error=mammalia_nt.%j.err      # STDERR
#SBATCH --partition=vm-small            # Partition/queue
#SBATCH --nodes=1                       # 1 node
#SBATCH --ntasks=1                      # 1 task
#SBATCH --cpus-per-task=16              # Threads for parallel
#SBATCH --time=05:00:00                 # hh:mm:ss
#SBATCH --account=OTH21013              # Your allocation

### 1) Load modules
module load biocontainers/0.1.0
module load hs-blastn
module load ncbi-blast-plus

export LC_ALL=C

### 2) Prepare BLASTDB directory
export BLASTDB=$WORK/Metabarcoding/blastdb
mkdir -p $BLASTDB && cd $BLASTDB

### 3) List and download NT volumes
rsync --list-only --no-motd rsync://ftp.ncbi.nlm.nih.gov/blast/db/nt_euk.*.tar.gz | awk '{print $NF}' > nt_euk.links.list             # list files :contentReference[oaicite:3]{index=3}

sed -i "s/'\$//" nt_euk.links.list

### 3.5) Download all volumes, resuming if interrupted
rsync --archive --partial --progress --compress --no-motd --files-from=nt_euk.links.list rsync://ftp.ncbi.nlm.nih.gov/blast/db/ .    # efficient batch download :contentReference[oaicite:9]{index=9}

### 4) Decompress 
for f in nt_euk.*.tar.gz; do
  tar -zxf "$f"
  rm *.tar.gz
done

### 5) Extract Mammalia sequences
blastdbcmd -db nt_euk -taxids 40674 -outfmt "%f" -out nt_euk_mammalia.fasta                                       # filter by taxid :contentReference[oaicite:6]{index=6}

### 6) Build Mammalia-only BLAST DB
makeblastdb -in nt_euk_mammalia.fasta -dbtype nucl -parse_seqids -out nt_euk_mammalia                                           # create custom DB :contentReference[oaicite:7]{index=7}

echo "Mammalia NT database ready in $BLASTDB/nt_euk_mammalia.*"
