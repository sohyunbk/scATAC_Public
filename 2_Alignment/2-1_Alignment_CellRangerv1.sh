#!/bin/bash
#SBATCH --job-name=Alignment        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=14             # Number of CPU cores per task
#SBATCH --mem=600gb                   # Job memory request #For normal fastq : 600gb
#SBATCH --time=120:00:04               # Time limit hrs:min:sec #For normal fastq : more than at least five days with 14 cpu 80 hours
#SBATCH --output=/scratch/sb14489/0.log/2_Mappingv1.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/2_Mappingv1.%j.err    # Standard error log
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)

Path=$1
Reference=/scratch/sb14489/0.Reference/$2
SampleName=$3
OutDirName=$4

module load CellRanger-ATAC/1.2.0
cd "$Path"/"$OutDirName"

cellranger-atac count \
   --id=$SampleName  \
   --reference=$Reference  \
   --fastqs=$Path/1.Rawdata/$SampleName  --localcores=14
