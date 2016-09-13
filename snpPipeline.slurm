#!/bin/sh
#SBATCH --time=160:00:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=24000	  # Maximum memory required per CPU (in megabyt$

#SBATCH --job-name=snps
#SBATCH --error=/work/hdzoo/rrculligan/cGAP_2.0/snps/snps.%J.err
#SBATCH --output=/work/hdzoo/rrculligan/cGAP_2.0/snps/snps.%J.out

#SBATCH -n 4

module load java/1.8
module load python/2.7
module load bwa/0.7;
module load compiler/gcc/4.8;
module load samtools/1.1;
module load bcftools/1.2-4-g1fedb8b;
module load tabix/0.2;
module load python/2.7;
module load novocraft;
module load blast-legacy/2.2.26

echo "Launching"
python snpPipe.py
echo "Done!"



