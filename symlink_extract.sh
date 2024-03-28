#!/bin/bash

# Source directory containing the files
SOURCE_DIR="/proj/uppmax2024-2-7/Genome_Analysis/4_Tean_Teh_2017/illumina_data"

# Target directory for the symbolic links
TARGET_DIR="/home/victoe/Genome_analysis/Data/RawData/PacBio"

# Pattern to match the desired files
PATTERN="*scaffold_10*.fastq.gz"

# Change to the source directory
cd "$SOURCE_DIR"

# Create symbolic links in the target directory for matching files
for file in $PATTERN; do
    ln -s "$(pwd)/$file" "$TARGET_DIR/"
done

echo "Symbolic links for '$PATTERN' have been created in '$TARGET_DIR'"

