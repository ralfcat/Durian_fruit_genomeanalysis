#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 14:00:00
#SBATCH -J eggnog_func_annotation
#SBATCH --mail-type=ALL
#SBATCH --mail-user=victor.englof.5352@student.uu.se
#SBATCH --output=eggnog_func_annotation.%j.out


module load bioinfo-tools
module load eggNOG-mapper/2.1.9
INPUT_DIR="/home/victoe/Genome_analysis/Data/Annotation/functional_annotation/output_proteins.fasta"
OUTPUT_DIR="/home/victoe/Genome_analysis/Data/Annotation/functional_annotation/eggnog_map"


# Run EggNOG-mapper
emapper.py -i ${INPUT_DIR}  --output ${OUTPUT_DIR} --override


echo "EggNOG-mapper run is complete. Check the output in ${OUTPUT_DIR}"
