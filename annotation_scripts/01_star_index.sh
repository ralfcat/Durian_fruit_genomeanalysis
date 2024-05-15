#!/bin/bash -l
#SBATCH -A uppmax2024-2-7            # Specify the project account for billing purposes.
#SBATCH -M snowy                      # Specifies the cluster name.
#SBATCH -p core                       # Specifies the partition (queue) your job will run in.
#SBATCH -n 8                          # Number of cores. Adjust based on your needs.
#SBATCH -t 04:00:00                   # Job runtime in HH:MM:SS format.
#SBATCH -J STAR_Index_masked      # Job name.
#SBATCH --mail-type=ALL               # Receive email notifications for job status.
#SBATCH --mail-user victor.englof.5352@student.uu.se  # Your email address.
#SBATCH --output=STAR_Index_masked.%j.out # Standard output and error log.

module load bioinfo-tools
module load star/2.7.11a

#set paths
GENOME_DIR=/home/victoe/Genome_analysis/Data/Annotation/star_index3 
GENOME_FASTA_FILES=/home/victoe/Genome_analysis/Data/Polishing/repeatmasker/pilon_polished.fasta.masked
READ_LENGTH=252  


#generate the STAR genome index
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir $GENOME_DIR \
     --genomeFastaFiles $GENOME_FASTA_FILES \



