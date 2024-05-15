#!/bin/bash

#source dir
SOURCE_DIR="/proj/uppmax2024-2-7/Genome_Analysis/4_Tean_Teh_2017/transcriptome/untrimmed"

#target directory for the symbolic links
TARGET_DIR="/home/victoe/Genome_analysis/Data/RawData/Illumina/RNAseq"

#pattern to match the desired files
PATTERN="*scaffold_10*.fastq.gz"


cd "$SOURCE_DIR"

for file in $PATTERN; do
    ln -s "$(pwd)/$file" "$TARGET_DIR/"
done

echo "Symbolic links for '$PATTERN' have been created in '$TARGET_DIR'"

