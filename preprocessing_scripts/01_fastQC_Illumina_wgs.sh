#!/bin/bash -l
#SBATCH -A uppmax2024-2-7        # Specify the project account for billing purposes.
#SBATCH -M snowy                 # Specifies the cluster name.
#SBATCH -p core                  # Specifies the partition (queue) your job will run in.
#SBATCH -n 2                     # Number of cores. Adjust if necessary, but do not exceed limits.
#SBATCH -t 00:20:00              # Job runtime in D-HH:MM:SS format (e.g., 02:00:00 for 2 hours).
#SBATCH -J FastQC_analysis       # Job name.
#SBATCH --mail-type=ALL          # Receive email notifications for job status.
#SBATCH --mail-user victor.englof.5352@student.uu.se # Your email address.
#SBATCH --output=%x.%j.out       # Standard output and error log.

module load bioinfo-tools
module load FastQC/0.11.9  

input_file="/home/victoe/Genome_analysis/Data/RawData/Illumina/WGS/SRR6058604_scaffold_10.2P.fastq.gz"
output_directory="/home/victoe/Genome_analysis/Data/PreprocessedData/QualityChecked"

#running fastqc
fastqc $input_file -o $output_directory -t 2

