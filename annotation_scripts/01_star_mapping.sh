#!/bin/bash -l
#SBATCH -A uppmax2024-2-7            # Specify the project account for billing purposes.
#SBATCH -M snowy                      # Specifies the cluster name.
#SBATCH -p core                       # Specifies the partition (queue) your job will run in.
#SBATCH -n 8                          # Number of cores. Adjust based on your needs.
#SBATCH -t 05:00:00                   # Job runtime in HH:MM:SS format.
#SBATCH -J STAR_Mapping               # Job name.
#SBATCH --mail-type=ALL               # Receive email notifications for job status.
#SBATCH --mail-user victor.englof.5352@student.uu.se  # Your email address.
#SBATCH --output=STAR_Mapping.%j.out  # Standard output and error log.

module load bioinfo-tools
module load star/2.7.11a

# Set paths
GENOME_DIR=/home/victoe/Genome_analysis/Data/Annotation/star_index
READS_DIR=/home/victoe/Genome_analysis/Data/PreprocessedData/TrimmedReads
OUTPUT_DIR=/home/victoe/Genome_analysis/Data/Annotation/star_mapping

read1_list=""
read2_list=""

for read1 in $READS_DIR/*_10.1.fastq.gz; do
	read1_list+="${read1},"
done

for read2 in $READS_DIR/*_10.2.fastq.gz; do
	read2_list+="${read2},"
done

#removes the last comma
read1_list=$(echo $read1_list | sed 's/,$//')
read2_list=$(echo $read2_list | sed 's/,$//')

STAR --runThreadN 8 \
     --genomeDir $GENOME_DIR \
     --readFilesIn $read1_list $read2_list \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_DIR/ \
     --outSAMtype BAM SortedByCoordinate \
