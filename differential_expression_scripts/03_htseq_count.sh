#!/bin/bash -l
#SBATCH -A uppmax2024-2-7            # Specify the project account for billing purposes.
#SBATCH -M snowy                      # Specifies the cluster name.
#SBATCH -p core                       # Specifies the partition (queue) your job will run in.
#SBATCH -n 2                          # Number of cores. Adjust based on your needs.
#SBATCH -t 05:00:00                   # Job runtime in HH:MM:SS format.
#SBATCH -J 03_htseq_count               # Job name.
#SBATCH --mail-type=ALL               # Receive email notifications for job status.
#SBATCH --mail-user victor.englof.5352@student.uu.se  # Your email address.
#SBATCH --output=htseq_cont.%j.out  # Standard output and error log.

module load bioinfo-tools
module load htseq/2.0.2

cd /home/victoe/Genome_analysis/Data/Annotation/star_mapping3



BAM_FILES=(./SRR*_scaffold_Aligned.sortedByCoord.out.bam)

# Count reads per gene using HTSeq
for BAM in "${BAM_FILES[@]}"; do
    echo "Processing $BAM"
    htseq-count -f bam -s no \
        $BAM /home/victoe/Genome_analysis/Data/Annotation/braker_output2/GeneMark-ET/genemark.gtf > $(basename $BAM .bam).counts2.txt
done


