#!/bin/bash -l
#SBATCH -A uppmax2024-2-7        # Specify the project account for billing purposes.
#SBATCH -M snowy                 # Specifies the cluster name.
#SBATCH -p core                  # Specifies the partition (queue) your job will run in.
#SBATCH -n 2                     # Number of cores. Adjust if necessary, but do not exceed limits.
#SBATCH -t 04:00:00              # Job runtime in D-HH:MM:SS format (e.g., 02:00:00 for 2 hours).
#SBATCH -J trimmomatic       # Job name.
#SBATCH --mail-type=ALL          # Receive email notifications for job status.
#SBATCH --mail-user victor.englof.5352@student.uu.se # Your email address.
#SBATCH --output=%x.%j.out       # Standard output and error log.

module load bioinfo-tools
module load trimmomatic/0.39

# Define your input/output directories
INPUT_DIR=/home/victoe/Genome_analysis/Data/RawData/Illumina/RNAseq
OUTPUT_DIR=/home/victoe/Genome_analysis/Data/RawData/Illumina/RNAseq/Trimmed_reads
ADAPTERS=/home/victoe/Genome_analysis/Code/preprocessing_scripts/TrueSeq2.fa

mkdir -p ${OUTPUT_DIR}

trimmomatic PE \
  -threads 2 \
  ${INPUT_DIR}/SRR6040095_scaffold_10.1.fastq.gz ${INPUT_DIR}/SRR6040095_scaffold_10.2.fastq.gz \
  ${OUTPUT_DIR}/sample1_1_paired.fastq.gz ${OUTPUT_DIR}/sample1_1_unpaired.fastq.gz \
  ${OUTPUT_DIR}/sample1_2_paired.fastq.gz ${OUTPUT_DIR}/sample1_2_unpaired.fastq.gz \
  ILLUMINACLIP:${ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
