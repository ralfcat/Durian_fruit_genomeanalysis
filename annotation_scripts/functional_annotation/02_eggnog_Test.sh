#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 14:00:00
#SBATCH -J eggnog_test_new
#SBATCH --mail-type=ALL
#SBATCH --mail-user=victor.englof.5352@student.uu.se
#SBATCH --output=eggnog_test_new.%j.out


module load bioinfo-tools
module load eggNOG-mapper/2.1.9
# module load gffread/0.12.7

ref_genome=/home/victoe/Genome_analysis/Data/Polishing/repeatmasker/pilon_polished.fasta.masked
genemark=/home/victoe/Genome_analysis/Data/Annotation/braker_output2/GeneMark-ET/genemark.gtf
output_dir=/home/victoe/Genome_analysis/Data/Annotation/eggmapper_new


# gffread $genemark -g $ref_genome -y ${output_dir}/predicted_proteins.faa
emapper.py -i ${output_dir}/predicted_proteins.faa -o ${output_dir}/emapper1 --override
