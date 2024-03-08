#!/bin/bash
#SBATCH --job-name=highmemtest        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=24            # Number of CPU cores per task
#SBATCH --mem=70gb                   # Job memory request
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/Check.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/Check.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail

module load  SAMtools/1.10-iccifort-2019.5.281
#cd /scratch/sb14489/3.scATAC/2.Maize_ear/2.Mapped_CellRangerv2/1_A619/outs
cd /scratch/sb14489/0.Reference/Maize_B73
samtools faidx Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa
