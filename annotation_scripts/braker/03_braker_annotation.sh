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

# Prepare AUGUSTUS and GeneMark environment
mkdir -p /home/victoe/Genome_analysis/Code/annotation_scripts/braker/augustus_config/species
chmod a+w -R /home/victoe/Genome_analysis/Code/annotation_scripts/braker/augustus_config/species/
cp -vf /sw/bioinfo/GeneMark/4.33-es/snowy/gm_key $HOME/.gm_key

# Load necessary modules
module load bioinfo-tools
module load braker/2.1.1_Perl5.24.1
module load augustus/3.2.3_Perl5.24.1
module load bamtools/2.5.1
module load blast/2.9.0+
module load GenomeThreader/1.7.0
module load samtools/1.8
module load GeneMark/4.33-es_Perl5.24.1

# Setup environment variables for AUGUSTUS and GeneMark
export AUGUSTUS_CONFIG_PATH=/home/victoe/Genome_analysis/Code/annotation_scripts/braker/augustus_config
export AUGUSTUS_BIN_PATH=/sw/bioinfo/augustus/3.4.0/snowy/bin
export AUGUSTUS_SCRIPTS_PATH=/sw/bioinfo/augustus/3.4.0/snowy/scripts
export GENEMARK_PATH=/sw/bioinfo/GeneMark/4.33-es/snowy

# Define paths for software and data
GENOME_PATH=/home/victoe/Genome_analysis/Data/Polishing/repeatmasker/pilon_polished.fasta.masked
OUT_DIR=/home/victoe/Genome_analysis/Data/Annotation/braker_output
mkdir -p $OUT_DIR

BAM_FILES=/home/victoe/Genome_analysis/Data/Annotation/star_mapping2

read1_list=""

for read1 in $BAM_FILES/*.bam; do
        read1_list+="${read1},"
done

read1_list=$(echo $read1_list | sed 's/,$//')

# Run BRAKER with all BAM files
braker.pl \
    --genome=$GENOME_PATH \
    --bam=$read1_list \
    --softmasking \
    --species="Durian" \
    --useexisting \
    --cores=8 \
    --workingdir=$OUT_DIR \
    --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH \
    --AUGUSTUS_BIN_PATH=$AUGUSTUS_BIN_PATH \
    --AUGUSTUS_SCRIPTS_PATH=$AUGUSTUS_SCRIPTS_PATH \
    --GENEMARK_PATH=$GENEMARK_PATH

