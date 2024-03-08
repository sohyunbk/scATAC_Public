#!/bin/bash
#SBATCH --job-name=MakeBW        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=12             # Number of CPU cores per task
#SBATCH --mem=200gb                   # Job memory request
#SBATCH --time=6:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/3_MakeBw_forGenomeBrowser.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/3_MakeBw_forGenomeBrowser.%j.err    # Standard error log

Path=$1
Reference=$2 #for fa.fai
SampleName=$3

cd $Path

## 1) Indexing
module load  SAMtools/1.10-iccifort-2019.5.281
samtools index -@ 24 "$SampleName"_Rmpcr.bam

## 2)Bigwig_forGenomeBrowser

module load Anaconda3/2020.02
source activate /home/sb14489/.conda/envs/ucsc
module load SAMtools/1.10-iccifort-2019.5.281
module load BEDTools/2.29.2-GCC-8.3.0

python /home/bth29393/jbscripts/file_to_bigwig_pe.py "$Reference"  "$SampleName"_Rmpcr.bam
bedtools bamtobed -i "$SampleName"_Rmpcr.bam > "$SampleName"_Rmpcr.bed
bedtools genomecov -i "$SampleName"_Rmpcr.bed -split -bg -g "$Reference"  > "$SampleName"_Rmpcr.bg
wigToBigWig "$SampleName"_Rmpcr.bg "$Reference"  "$SampleName"_Rmpcr.bw
