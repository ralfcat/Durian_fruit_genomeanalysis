#!/bin/bash -l
#SBATCH -A uppmax2024-2-7        # Specify the project account for billing purposes.
#SBATCH -M snowy                 # Specifies the cluster name.
#SBATCH -p core                  # Specifies the partition (queue) your job will run in.
#SBATCH -n 4                     # Adjust the number of cores according to your needs and Flye's requirements.
#SBATCH -t 6:00:00            # Adjust job runtime as needed in D-HH:MM:SS format (e.g., 1-00:00:00 for 1 day).
#SBATCH -J Flye_Assembly         # Job name.
#SBATCH --mail-type=ALL          # Receive email notifications for job status.
#SBATCH --mail-user victor.englof.5352@student.uu.se  # Your email address.
#SBATCH --output=Flye_Assembly.%j.out  # Standard output and error log.


module load bioinfo-tools
module load Flye/2.9.1          

#run the assembly
flye --pacbio-raw /home/victoe/Genome_analysis/Data/RawData/PacBio/SRR6037732_scaffold_10.fq.gz \
     --out-dir /home/victoe/Genome_analysis/Data/Assembly/PacBioAssembly \
     --threads 4                 # Match this to the number of cores requested with -n.
