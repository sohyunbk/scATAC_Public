#!/bin/bash
#SBATCH --job-name=Alignment        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=14             # Number of CPU cores per task
#SBATCH --mem=600gb                   # Job memory request #For normal fastq : 600gb
#SBATCH --time=30:00:04               # Time limit hrs:min:sec #For normal fastq : 80 hours
#SBATCH --output=/scratch/sb14489/0.log/2_Mapping.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/2_Mapping.%j.err    # Standard error log
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)

Path=$1
Reference=/scratch/sb14489/0.Reference/$2
SampleName=$3
OutDirName=$4

module load CellRanger-ATAC/2.0.0
cd "$Path"/"$OutDirName"

cellranger-atac count \
   --id=$SampleName  \
   --reference=$Reference  \
   --fastqs=$Path/1.Rawdata/$SampleName  --localcores=14

## They search the subdirectory for rawdata too.
