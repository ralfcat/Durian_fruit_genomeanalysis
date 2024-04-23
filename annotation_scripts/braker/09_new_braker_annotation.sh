#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 12:00:00
#SBATCH -J BRAKER_annotation
#SBATCH --mail-type=ALL
#SBATCH --mail-user=victor.englof.5352@student.uu.se
#SBATCH --output=BRAKER_annotation.%j.out

# Load required modules
module load bioinfo-tools
module load braker/2.1.1_Perl5.24.1
module load augustus/3.2.3_Perl5.24.1
module load bamtools/2.5.1
module load blast/2.9.0+
module load GenomeThreader/1.7.0
module load samtools/1.8
module load GeneMark/4.33-es_Perl5.24.1

OUT_DIR=/home/victoe/Genome_analysis/Data/Annotation/braker_output2

# Set environment variables
export AUGUSTUS_CONFIG_PATH=/home/victoe/Genome_analysis/Data/Annotation/braker_output2/augustus_config
export AUGUSTUS_BIN_PATH=/sw/bioinfo/augustus/3.4.0/snowy/bin
export AUGUSTUS_SCRIPTS_PATH=/sw/bioinfo/augustus/3.4.0/snowy/scripts
export GENEMARK_PATH=/sw/bioinfo/GeneMark/4.33-es/snowy

# Update AUGUSTUS configuration permissions
mkdir -p ${OUT_DIR}/augustus_config/species
chmod a+w -R ${OUT_DIR}/augustus_config/species

# Copy GeneMark key
cp -vf /sw/bioinfo/GeneMark/4.33-es/snowy/gm_key $HOME/.gm_key

# Copy AUGUSTUS configuration file
source $AUGUSTUS_CONFIG_COPY
echo "AUGUSTUS_CONFIG_COPY is set to: $AUGUSTUS_CONFIG_COPY"


# Replace the following paths with the correct paths to your input files
GENOME="/home/victoe/Genome_analysis/Data/Polishing/repeatmasker/pilon_polished.fasta.masked"
BAM=$(ls /home/victoe/Genome_analysis/Data/Annotation/star_mapping3/*_scaffold_Aligned.sortedByCoord.out.bam | paste -sd, -)


# Running BRAKER
braker.pl \
   --genome=${GENOME} \
   --bam=${BAM} \
   --species=your_species \
   --softmasking \
   --useexisting \
   --cores=8	\
   --workingdir=${OUT_DIR}
