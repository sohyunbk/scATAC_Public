#!/bin/bash
#SBATCH --job-name=Picard        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=4             # Number of CPU cores per task
#SBATCH --mem=400gb                   # Job memory request
#SBATCH --time=12:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/4_Picards_Class_v16.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/4_Picards_Class_v16.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail

## it can have memory error if it's too low Runtime.totalMemory()=2058354688
Path=$1
SampleName=$2

cd "$Path"/3.SortedBam
#cd "$Path"/3.SortedBam_Class
conda activate /home/sb14489/.conda/envs/r_env

module load picard/2.16.0-Java-1.8.0_144
#module load picard/2.21.6-Java-11
#module load picard/2.26.10-Java-13

#java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
#    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100 REMOVE_DUPLICATES=false METRICS_FILE= ./"$SampleName"_dups_Markingpcr.txt \
#    I= ./"$SampleName"_Sorted.bam \
#    O= ./"$SampleName"_Markingpcr.bam \
#    BARCODE_TAG=CB \
#    ASSUME_SORT_ORDER=coordinate\
#    USE_JDK_DEFLATER=true USE_JDK_INFLATER=true TMP_DIR="/scratch/sb14489/0.log/"

java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100 REMOVE_DUPLICATES=true METRICS_FILE= ./"$SampleName"_dups.txt \
        I= ./"$SampleName"_Sorted.bam \
        O= ./"$SampleName"_Rmpcr.bam \
        BARCODE_TAG=CB \
        ASSUME_SORT_ORDER=coordinate\
        USE_JDK_DEFLATER=true USE_JDK_INFLATER=true TMP_DIR="/scratch/sb14489/0.log/"
