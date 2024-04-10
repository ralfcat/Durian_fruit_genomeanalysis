#!/bin/bash -l
#SBATCH -A uppmax2024-2-7          # Specify the project account for billing purposes.
#SBATCH -M snowy                    # Specifies the cluster name.
#SBATCH -p core                     # Specifies the partition (queue) your job will run in.
#SBATCH -n 2                       # Number of cores. Adjust based on your needs.
#SBATCH -t 02:00:00               # Job runtime in D-HH:MM:SS format (e.g., 1-00:00:00 for 1 day).
#SBATCH -J Pilon_polishing          # Job name.
#SBATCH --mail-type=ALL             # Receive email notifications for job status.
#SBATCH --mail-user victor.englof.5352@student.uu.se  # Your email address.
#SBATCH --output=Pilon_polishing.%j.out # Standard output and error log.

module load bioinfo-tools
module load Pilon/1.24 

# Set paths
ASSEMBLY_PATH=/home/victoe/Genome_analysis/Data/Assembly/PacBioAssembly/assembly.fasta
BAM_DIR=/home/victoe/Genome_analysis/IlluminaAssembly/BWA_alignment/SRR6058604_scaffold_10.sorted.bam
OUT_DIR=/home/victoe/Genome_analysis/Data/Polishing/Pilon_polished
mkdir -p $OUT_DIR



java -Xmx14G -jar $PILON_HOME/pilon.jar \
    --genome $ASSEMBLY_PATH \
    --frags $BAM_DIR \
    --output $OUT_DIR/pilon_polished \
    --threads 2 \
    --changes



