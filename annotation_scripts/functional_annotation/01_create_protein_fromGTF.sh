#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 1:00:00
#SBATCH -J create_fasta
#SBATCH --mail-type=ALL
#SBATCH --mail-user=victor.englof.5352@student.uu.se
#SBATCH --output=create_fasta.%j.out

module load bioinfo-tools
module load cufflinks/2.2.1
module load emboss/6.6.0

GENOME_PATH="/home/victoe/Genome_analysis/Data/Polishing/repeatmasker/pilon_polished.fasta.masked"
GTF_PATH="/home/victoe/Genome_analysis/Data/Annotation/braker_output2/GeneMark-ET/genemark.gtf"
OUTPUT_DIR="/home/victoe/Genome_analysis/Data/Annotation/functional_annotation"

# Extract CDS sequences from GTF using gffread
gffread -x ${OUTPUT_DIR}/output_cds.fasta -g ${GENOME_PATH} ${GTF_PATH}

# Translate CDS to protein sequences using transeq
transeq ${OUTPUT_DIR}/output_cds.fasta ${OUTPUT_DIR}/output_proteins.fasta -frame=1 -clean

echo "Extraction and translation complete. Results are in ${OUTPUT_DIR}"

