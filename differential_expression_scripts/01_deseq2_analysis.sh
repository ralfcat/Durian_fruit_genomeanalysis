#!/bin/bash -l
#SBATCH -A uppmax2024-2-7            # Specify the project account for billing purposes.
#SBATCH -M snowy                      # Specifies the cluster name.
#SBATCH -p core                       # Specifies the partition (queue) your job will run in.
#SBATCH -n 2                          # Number of cores. Adjust based on your needs.
#SBATCH -t 05:00:00                   # Job runtime in HH:MM:SS format.
#SBATCH -J 01_deseq2                  # Job name.
#SBATCH --mail-type=ALL               # Receive email notifications for job status.
#SBATCH --mail-user victor.englof.5352@student.uu.se  # Your email address.
#SBATCH --output=deseq2_cont.%j.out  # Standard output and error log.

module load bioinfo-tools
module load R_packages/4.3.1

cd /home/victoe/Genome_analysis/Code/differential_expression_scripts

Rscript 02_deseq2.R

