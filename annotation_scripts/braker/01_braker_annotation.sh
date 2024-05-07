#!/bin/bash -l
#SBATCH -A uppmax2024-2-7            # Specify the project account for billing purposes.
#SBATCH -M snowy                      # Specifies the cluster name.
#SBATCH -p core                       # Specifies the partition (queue) your job will run in.
#SBATCH -n 8                          # Number of cores. Adjust based on your needs.
#SBATCH -t 18:00:00                 # Job runtime in D-HH:MM:SS format.
#SBATCH -J BRAKER_annotation          # Job name.
#SBATCH --mail-type=ALL               # Receive email notifications for job status.
#SBATCH --mail-user victor.englof.5352@student.uu.se  # Your email address.
#SBATCH --output=BRAKER_annotation.%j.out # Standard output and error log.

# Load necessary modules
module load bioinfo-tools
module load braker/2.1.1_Perl5.24.1
module load augustus/3.2.3_Perl5.24.1
module load bamtools/2.5.1
module load blast/2.9.0+
module load GenomeThreader/1.7.0
module load samtools/1.8
module load GeneMark/4.33-es_Perl5.24.1

export AUGUSTUS_CONFIG_PATH=/home/victoe/augustus_config
source $AUGUSTUS_CONFIG_COPY

cp -vf /sw/bioinfo/GeneMark/4.33-es/snowy/gm_key $HOME/.gm_key

chmod a+w -R /home/victoe/augustus_config/species/

#define paths for software
export AUGUSTUS_BIN_PATH=/sw/bioinfo/augustus/3.4.0/snowy/bin
export AUGUSTUS_SCRIPTS_PATH=/sw/bioinfo/augustus/3.4.0/snowy/scripts
export GENEMARK_PATH=/sw/bioinfo/GeneMark/4.33-es/snowy

#set paths for your genome and RNA-Seq BAM files
GENOME_PATH=/home/victoe/Genome_analysis/Data/Polishing/repeatmasker/pilon_polished.fasta.masked
BAM_FILE=/path/to/your/rnaseq.bam
OUT_DIR=/path/to/your/braker_output
mkdir -p $OUT_DIR

# Run BRAKER
braker.pl \
    --genome=$GENOME_PATH \
    --bam=$BAM_FILE \
    --softmasking \
    --useexisting \
    --cores=8 \
    --species=your_species_name \
    --workingdir=$OUT_DIR \
    --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH \
    --AUGUSTUS_BIN_PATH=$AUGUSTUS_BIN_PATH \
    --AUGUSTUS_SCRIPTS_PATH=$AUGUSTUS_SCRIPTS_PATH \
    --GENEMARK_PATH=$GENEMARK_PATH

echo "BRAKER annotation completed."

