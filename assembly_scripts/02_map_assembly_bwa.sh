#!/bin/bash -l
#SBATCH -A uppmax2024-2-7            
#SBATCH -M snowy                      
#SBATCH -p core                       
#SBATCH -n 2                          
#SBATCH -t 02:00:00                 
#SBATCH -J BWA_alignment              
#SBATCH --mail-type=ALL               
#SBATCH --mail-user victor.englof.5352@student.uu.se  
#SBATCH --output=BWA_alignment.%j.out 


module load bioinfo-tools
module load bwa/0.7.17                
module load samtools/1.10             

# path to assembly
ASSEMBLY_PATH=/home/victoe/Genome_analysis/Data/Assembly/PacBioAssembly/assembly.fasta

# Check if the assembly is indexed; if not, index it
if [ ! -f $ASSEMBLY_PATH.bwt ]; then
    bwa index $ASSEMBLY_PATH
fi

#path to read files
READ1=/home/victoe/Genome_analysis/Data/RawData/Illumina/WGS/SRR6058604_scaffold_10.1P.fastq.gz
READ2=/home/victoe/Genome_analysis/Data/RawData/Illumina/WGS/SRR6058604_scaffold_10.2P.fastq.gz

#output dir
OUT_DIR=/home/victoe/Genome_analysis/Data/Assembly/IlluminaAssembly/BWA_alignment
mkdir -p $OUT_DIR

# output file name
SAMPLE_NAME="SRR6058604_scaffold_10"
BAM_FILE=$OUT_DIR/${SAMPLE_NAME}.sorted.bam

# Perform alignment, convert to BAM, sort, and index
bwa mem -t 2 $ASSEMBLY_PATH $READ1 $READ2 | \
samtools view -@ 2 -bu - | \
samtools sort -@ 2 -o $BAM_FILE

samtools index $BAM_FILE




