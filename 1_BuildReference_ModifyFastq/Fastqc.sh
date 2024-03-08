#!/bin/bash
#SBATCH --job-name=LowMemory        # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=100gb                     # Job memory request
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/3.scATAC_flo/100.log/fastqc.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/3.scATAC_flo/100.log/fastqc.%j.err    # Standard error log
#SBATCH --array=0-2                   # Array range

#List=(1_A619_2  2_rel2_2  3_bif3_2  4_relk1_2)
List=(Tassel_ba1 Tassel_bif1 Seedling_IAA2)

Path=/scratch/sb14489/3.scATAC_flo

module load FastQC/0.11.9-Java-11

cd $Path/1.Rawdata/"${List[SLURM_ARRAY_TASK_ID]}"

#fastqc -t 32 *

R1File="*R1*.zip"
echo $R1File 
module load UnZip/6.0-GCCcore-10.3.0
unzip $R1File
