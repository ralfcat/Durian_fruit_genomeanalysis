#!/bin/bash -l
#SBATCH -A uppmax2024-2-7            
#SBATCH -M snowy                      
#SBATCH -p core                       
#SBATCH -n 8                          
#SBATCH -t 18:00:00                 
#SBATCH -J BRAKER_annotation          
#SBATCH --mail-type=ALL               
#SBATCH --mail-user victor.englof.5352@student.uu.se  
#SBATCH --output=BRAKER_annotation.%j.out 

# Load necessary modules
module load bioinfo-tools
module load braker/2.1.1_Perl5.24.1
module load augustus/3.2.3_Perl5.24.1
module load bamtools/2.5.1
module load blast/2.9.0+
module load GenomeThreader/1.7.0
module load samtools/1.8
module load GeneMark/4.33-es_Perl5.24.1

# Ensure the gm_key is copied and accessible
cp -vf /sw/bioinfo/GeneMark/4.33-es/snowy/gm_key $HOME/.gm_key
chmod 600 $HOME/.gm_key

# Ensure AUGUSTUS config directory is local and writable
mkdir -p $HOME/augustus_config
export AUGUSTUS_CONFIG_PATH=$HOME/augustus_config

# Ensure permissions are set correctly for AUGUSTUS to write into the config directory
chmod -R 775 $AUGUSTUS_CONFIG_PATH

# Ensure paths are properly defined for executables
export AUGUSTUS_BIN_PATH=/sw/bioinfo/augustus/3.4.0/snowy/bin
export AUGUSTUS_SCRIPTS_PATH=/sw/bioinfo/augustus/3.4.0/snowy/scripts
export GENEMARK_PATH=/sw/bioinfo/GeneMark/4.33-es/snowy

# Define paths for software and data
GENOME_PATH=/home/victoe/Genome_analysis/Data/Polishing/repeatmasker/pilon_polished.fasta.masked
OUT_DIR=/home/victoe/Genome_analysis/Data/Annotation/braker_output
mkdir -p $OUT_DIR

# Combine all BAM files into one string with commas
BAM_FILES=$(ls /home/victoe/Genome_analysis/Data/Annotation/star_mapping2/*.bam | tr '\n' ',' | sed 's/,$//')

# Run BRAKER with all BAM files
braker.pl \
    --genome=$GENOME_PATH \
    --bam=$BAM_FILES \
    --softmasking \
    --species="Durian" \
    --useexisting \
    --cores=8 \
    --workingdir=$OUT_DIR \
    --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH \
    --AUGUSTUS_BIN_PATH=$AUGUSTUS_BIN_PATH \
    --AUGUSTUS_SCRIPTS_PATH=$AUGUSTUS_SCRIPTS_PATH \
    --GENEMARK_PATH=$GENEMARK_PATH

