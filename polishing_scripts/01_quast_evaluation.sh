#!/bin/bash -l
#SBATCH -A uppmax2024-2-7          # Specify the project account for billing purposes.
#SBATCH -M snowy                    # Specifies the cluster name.
#SBATCH -p core                     # Specifies the partition (queue) your job will run in.
#SBATCH -n 1                       # Number of cores. Adjust based on your needs.
#SBATCH -t 02:00:00               # Job runtime in D-HH:MM:SS format (e.g., 1-00:00:00 for 1 day).
#SBATCH -J quast_eval          # Job name.
#SBATCH --mail-type=ALL             # Receive email notifications for job status.
#SBATCH --mail-user victor.englof.5352@student.uu.se  # Your email address.
#SBATCH --output=quast_eval.%j.out # Standard output and error log.

module load bioinfo-tools
module load quast/5.0.2

POLISHED_ASSEMBLY_PATH=/home/victoe/Genome_analysis/Data/Polishing/Pilon_polished/pilon_polished.fasta
OUT_DIR=/home/victoe/Genome_analysis/Data/Polishing/quast_eval
FLYE_ASSEMBLY_PATH=/home/victoe/Genome_analysis/Data/Assembly/PacBioAssembly/assembly.fasta

#running quast
quast.py $POLISHED_ASSEMBLY_PATH  -o $OUT_DIR

