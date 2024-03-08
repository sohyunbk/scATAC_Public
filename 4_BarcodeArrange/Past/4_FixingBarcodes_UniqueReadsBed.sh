#!/bin/bash
#SBATCH --job-name=FixingBarcode        # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=80gb                   # Job memory request
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/4_FixingBarcodes.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/4_FixingBarcodes.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail

Path=$1 #/scratch/sb14489/3.scATAC_flo
ConfigFile=$2
SampleName=$3

cd "$Path"/3.SortedBam
mkdir -p "$Path"/4.Bam_FixingBarcode

module load Anaconda3/2020.02
source activate /home/sb14489/.conda/envs/ucsc
module load  SAMtools/1.10-iccifort-2019.5.281

#FixingBarcode
python /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/4_BarcodeArrange/fixBC_v2.py -10x_config \
"$ConfigFile"  -BAM "$SampleName"_Rmpcr.bam  \
-exp_name "$SampleName"  \
-tissue Temp \
-output ../4.Bam_FixingBarcode/"$SampleName".count_file_metrics | samtools view -bh - > ../4.Bam_FixingBarcode/"$SampleName"_BarcodeFixed.bam

#~/.conda/env/bin/python3.7 /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/4_BarcodeArrange/fixBC_v2.py -10x_config \
#/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai  -BAM 4_relk1_2_Rmpcr.bam \
#-exp_name 4_relk1_2  \
#-tissue Temp \
#-output ../4.Bam_FixingBarcode/4_relk1_2.count_file_metrics | samtools view -bh - > ../4.Bam_FixingBarcode/4_relk1_2_BarcodeFixed.bam



#UniqueReadsBed
perl /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/4_BarcodeArrange/makeTn5bed.pl ../4.Bam_FixingBarcode/"$SampleName"_BarcodeFixed.bam  | sort -k1,1 -k2,2n - > ../4.Bam_FixingBarcode/"$SampleName"_Unique.bed

#TAGCCCTCACATTGCA
#ATCGGGAGTGTAACGT
#TGCAATGTGAGGGCTA
