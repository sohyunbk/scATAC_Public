#!/bin/bash
#SBATCH --job-name=FilteringREads        # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=12            # Number of CPU cores per task
#SBATCH --mem=90gb                   # Job memory request
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/3_RemoveMultiMappingReads_2ndRe.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/3_RemoveMultiMappingReads_2ndRe.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail

## It's for the output of cellranger v2
Path=$1
MappedDir=$2
SampleName=$3

module load  SAMtools/1.10-iccifort-2019.5.281
cd "$Path"
mkdir -p 3.SortedBam

##########################################################################################
## Filtering
##########################################################################################
samtools view -@ 12 -h -f 3 -q 10 "$Path"/"$MappedDir"/"$SampleName"/outs/possorted_bam.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -@ 12 -bS - > "$Path"/3.SortedBam/"$SampleName"_Sorted.bam
#-h, --with-header          Include header in SAM output
#-f, --require-flags FLAG   ...have all of the FLAGs present
# -q, --min-MQ INT           ...have mapping quality >= INT
# -v -e : exclude the lines : XA:Z: SA:Z:(rname ,pos ,strand ,CIGAR ,mapQ ,NM ;)+ Other canonical alignments in a chimeric alignment, formatted as a semicolon-delimited list. Each element in the list represents a part of the chimeric alignment. Conventionally, at a supplementary line, the first element points to the primary line. Strand is either ‘+’ or ‘-’, indicating forward/reverse strand, corresponding to FLAG bit 0x10. Pos is a 1-based coordinate.

# XA: Alternative hits; format: (chr,pos,CIGAR,NM;)*
# SA: Z, not sure what it is but, it almost always coincides with the 256 flag = not primary alignment

## Edited Version
#samtools view -@ 32 -h -f 3 -q 20 "$Path"/"$MappedDir"/"$SampleName"/outs/possorted_bam.bam | samtools view -@ 12 -bS - > "$Path"/3.SortedBam/"$SampleName"_Sorted.bam

# E.g.
# sbatch 3-1_RemoveMultiMappingReads.sh /scratch/sb14489/3.scATAC/2.Maize_ear 2.Mapped_CellRangerv2 1_A619
# Check the fragmentsize
#awk '$9 > 2000 || $9 < -2000' 3_bif3_BarcodeFixed.sam 
##########################################################################################
## Check and change Header for PICARD and Cellrangerv2
##########################################################################################
samtools view -H "$Path"/3.SortedBam/"$SampleName"_Sorted.bam > "$Path"/3.SortedBam/"$SampleName"_Sorted_Header.sam
sed -i '1d' "$Path"/3.SortedBam/"$SampleName"_Sorted_Header.sam
sed -i '1i @HD\tVN:1.5\tSO:coordinate' "$Path"/3.SortedBam/"$SampleName"_Sorted_Header.sam
samtools reheader "$Path"/3.SortedBam/"$SampleName"_Sorted_Header.sam  "$Path"/3.SortedBam/"$SampleName"_Sorted.bam > "$Path"/3.SortedBam/"$SampleName"_Sorted_HF.bam
