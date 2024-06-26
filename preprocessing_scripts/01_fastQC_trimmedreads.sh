#!/bin/bash -l
#SBATCH -A uppmax2024-2-7        # Specify the project account for billing purposes.
#SBATCH -M snowy                 # Specifies the cluster name.
#SBATCH -p core                  # Specifies the partition (queue) your job will run in.
#SBATCH -n 2                     # Number of cores. Adjust if necessary, but do not exceed limits.
#SBATCH -t 00:20:00              # Job runtime in D-HH:MM:SS format (e.g., 02:00:00 for 2 hours).
#SBATCH -J FastQC_trimmed       # Job name.
#SBATCH --mail-type=ALL          # Receive email notifications for job status.
#SBATCH --mail-user victor.englof.5352@student.uu.se # Your email address.
#SBATCH --output=%x.%j.out       # Standard output and error log.

# Load required modules
module load bioinfo-tools
module load FastQC/0.11.9  # Assuming FastQC is the correct module name. Adjust as necessary.

#Run FastQC on your desired file with specified number of threads
# Replace 'yourfile.fastq.gz' with the path to your actual file and 'output_directory/' with your desired output directory.
input_file="/home/victoe/Genome_analysis/Data/RawData/Illumina/RNAseq/Trimmed_reads/sample1_1_paired.fastq.gz"
output_directory="/home/victoe/Genome_analysis/Data/RawData/Illumina/RNAseq/Trimmed_reads"

fastqc $input_file -o $output_directory -t 2

