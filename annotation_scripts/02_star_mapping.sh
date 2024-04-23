#!/bin/bash -l
#SBATCH -A uppmax2024-2-7            # Specify the project account for billing purposes.
#SBATCH -M snowy                      # Specifies the cluster name.
#SBATCH -p core                       # Specifies the partition (queue) your job will run in.
#SBATCH -n 8                          # Number of cores. Adjust based on your needs.
#SBATCH -t 05:00:00                   # Job runtime in HH:MM:SS format.
#SBATCH -J 01_STAR_Mapping               # Job name.
#SBATCH --mail-type=ALL               # Receive email notifications for job status.
#SBATCH --mail-user victor.englof.5352@student.uu.se  # Your email address.
#SBATCH --output=02_STAR_Mapping.%j.out  # Standard output and error log.

module load bioinfo-tools
module load star/2.7.11a

GENOME_DIR=/home/victoe/Genome_analysis/Data/Annotation/star_index3
READS_DIR=/home/victoe/Genome_analysis/Data/PreprocessedData/TrimmedReads
OUTPUT_DIR=/home/victoe/Genome_analysis/Data/Annotation/star_mapping3

mkdir -p $OUTPUT_DIR

for read1 in $READS_DIR/*_10.1.fastq.gz; do
    read2="${read1/_10.1.fastq.gz/_10.2.fastq.gz}"  
    sample_name=$(basename $read1 _10.1.fastq.gz)  
    
    STAR --runThreadN 8 \
         --genomeDir $GENOME_DIR \
         --readFilesIn $read1 $read2 \
         --readFilesCommand zcat \
         --outFileNamePrefix $OUTPUT_DIR/${sample_name}_ \
         --outSAMtype BAM SortedByCoordinate
done

