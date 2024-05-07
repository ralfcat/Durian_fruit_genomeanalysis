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

# Prepare GeneMark by copying the key and setting permissions
cp -vf /sw/bioinfo/GeneMark/4.33-es/snowy/gm_key $HOME/.gm_key
chmod 600 $HOME/.gm_key

# Load necessary modules
module load bioinfo-tools
module load braker/2.1.1_Perl5.24.1
module load augustus/3.2.3_Perl5.24.1
module load bamtools/2.5.1
module load blast/2.9.0+
module load GenomeThreader/1.7.0
module load samtools/1.8
module load GeneMark/4.33-es_Perl5.24.1

# Define paths for software and data
GENOME_PATH=/home/victoe/Genome_analysis/Data/Polishing/repeatmasker/pilon_polished.fasta.masked
OUT_DIR=/home/victoe/Genome_analysis/Data/Annotation/braker_output
mkdir -p $OUT_DIR

# Directory containing the BAM files
BAM_DIR=/home/victoe/Genome_analysis/Data/Annotation/star_mapping2

# Check BAM file accessibility and permissions
echo "Validating BAM files..."
for bam in $BAM_DIR/*.bam; do
    if [ -s "$bam" ] && [ -r "$bam" ]; then  # Check if file is not empty and is readable
        echo "Valid BAM: $bam"
    else
        echo "Invalid or inaccessible BAM file skipped: $bam"
    fi
done

# Creating a comma-separated list of BAM files (ensuring readability and correctness)
BAM_FILES=$(find $BAM_DIR -name "*.bam" -size +1c -print0 | tr '\0' ',' | sed 's/,$//')

if [ -z "$BAM_FILES" ]; then
    echo "No valid BAM files found. Exiting."
    exit 1
fi

echo "Using BAM files: $BAM_FILES"

# Run BRAKER with all valid BAM files
braker.pl \
    --genome=$GENOME_PATH \
    --bam=$BAM_FILES \
    --softmasking \
    --species="Durian" \
    --useexisting \
    --cores=8 \
    --workingdir=$OUT_DIR \
    --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH \
    --AUGUSTUS_BIN_PATH=/sw/bioinfo/augustus/3.4.0/snowy/bin \
    --AUGUSTUS_SCRIPTS_PATH=/sw/bioinfo/augustus/3.4.0/snowy/scripts \
    --GENEMARK_PATH=$GENEMARK_PATH

