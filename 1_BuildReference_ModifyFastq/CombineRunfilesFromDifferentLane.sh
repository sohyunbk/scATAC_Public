#!/bin/bash
#SBATCH --job-name=CombineFiles       # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=80gb                     # Job memory request
#SBATCH --time=6:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/3.scATAC_flo/100.log/1_CombineFiles.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/3.scATAC_flo/100.log/1_CombineFiles.%j.err    # Standard error log

Path=/scratch/sb14489/3.scATAC_flo/1.Rawdata/

cd $Path

#mkdir 3_bif3_2
#cat ./3_bif3_2ndRun_2Re/3_bif3_2_S6_L003_I1_001.fastq.gz ./3_bif3_3rdRun_2Re/3_bif3_3_S6_L002_I1_001.fastq.gz > ./3_bif3_2/3_bif3_2_S6_L002_I1_001.fastq.gz
#cat ./3_bif3_2ndRun_2Re/3_bif3_2_S6_L003_R1_001.fastq.gz ./3_bif3_3rdRun_2Re/3_bif3_3_S6_L002_R1_001.fastq.gz > ./3_bif3_2/3_bif3_2_S6_L002_R1_001.fastq.gz
#cat ./3_bif3_2ndRun_2Re/3_bif3_2_S6_L003_R2_001.fastq.gz ./3_bif3_3rdRun_2Re/3_bif3_3_S6_L002_R2_001.fastq.gz > ./3_bif3_2/3_bif3_2_S6_L002_R2_001.fastq.gz
#cat ./3_bif3_2ndRun_2Re/3_bif3_2_S6_L003_R3_001.fastq.gz ./3_bif3_3rdRun_2Re/3_bif3_3_S6_L002_R3_001.fastq.gz > ./3_bif3_2/3_bif3_2_S6_L002_R3_001.fastq.gz

#mkdir 4_relk1_2
#cat ./4_relk1_2ndRun_2Re/4_relk1_2_S6_L001_I1_001.fastq.gz ./4_relk1_3rdRun_2Re/4_relk1_3_S5_L003_I1_001.fastq.gz > ./4_relk1_2/4_relk1_2_S5_L003_I1_001.fastq.gz
#cat ./4_relk1_2ndRun_2Re/4_relk1_2_S6_L001_R1_001.fastq.gz  ./4_relk1_3rdRun_2Re/4_relk1_3_S5_L003_R1_001.fastq.gz > ./4_relk1_2/4_relk1_2_S5_L003_R1_001.fastq.gz
#cat ./4_relk1_2ndRun_2Re/4_relk1_2_S6_L001_R2_001.fastq.gz  ./4_relk1_3rdRun_2Re/4_relk1_3_S5_L003_R2_001.fastq.gz > ./4_relk1_2/4_relk1_2_S5_L003_R2_001.fastq.gz
#cat ./4_relk1_2ndRun_2Re/4_relk1_2_S6_L001_R3_001.fastq.gz  ./4_relk1_3rdRun_2Re/4_relk1_3_S5_L003_R3_001.fastq.gz > ./4_relk1_2/4_relk1_2_S5_L003_R3_001.fastq.gz

cat ./Tassel_WT/Tassel_S1_L001_I1_001.fastq.gz ./Tassel_WT/Tassel_S1_L002_I1_001.fastq.gz ./Tassel_WT/Tassel_S1_L003_I1_001.fastq.gz ./Tassel_WT/Tassel_S1_L004_I1_001.fastq.gz > ./Tassel_WT/Tassel_S1_L111_I1_001.fastq.gz
cat ./Tassel_WT/Tassel_S1_L001_R1_001.fastq.gz ./Tassel_WT/Tassel_S1_L002_R1_001.fastq.gz ./Tassel_WT/Tassel_S1_L003_R1_003.fastq.gz ./Tassel_WT/Tassel_S1_L004_R1_001.fastq.gz > ./Tassel_WT/Tassel_S1_L111_R1_001.fastq.gz
cat ./Tassel_WT/Tassel_S1_L001_R2_001.fastq.gz  ./Tassel_WT/Tassel_S1_L002_R2_001.fastq.gz  ./Tassel_WT/Tassel_S1_L003_R2_001.fastq.gz  ./Tassel_WT/Tassel_S1_L004_R2_001.fastq.gz >  ./Tassel_WT/Tassel_S1_L111_R2_001.fastq.gz
cat ./Tassel_WT/Tassel_S1_L001_R3_001.fastq.gz  ./Tassel_WT/Tassel_S1_L002_R3_001.fastq.gz  ./Tassel_WT/Tassel_S1_L003_R3_001.fastq.gz  ./Tassel_WT/Tassel_S1_L004_R3_001.fastq.gz >  ./Tassel_WT/Tassel_S1_L111_R3_001.fastq.gz
