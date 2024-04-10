#!/bin/bash -l
#SBATCH -A uppmax2024-2-7            
#SBATCH -M snowy                      
#SBATCH -p core                       
#SBATCH -n 4                          
#SBATCH -t 02:00:00                 
#SBATCH -J repeatmasker_polishing              
#SBATCH --mail-type=ALL               
#SBATCH --mail-user victor.englof.5352@student.uu.se  
#SBATCH --output=01_.%j.out

module load bioinfo-tools
module load RepeatMasker/4.1.5

assembly_path=/home/victoe/Genome_analysis/Data/Polishing/Pilon_polished/pilon_polished.fasta
out_path=/home/victoe/Genome_analysis/Data/Polishing/repeatmasker
 
RepeatMasker -pa 4 -dir $out_path $assembly_path


