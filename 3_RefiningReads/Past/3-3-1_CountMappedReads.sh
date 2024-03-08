#!/bin/bash
#SBATCH --job-name=MappingRate        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=12             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/3.scATAC_flo/100.log/3-3_CountMappedReads_%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/3.scATAC_flo/100.log/3-3_CountMappedReads_%j.err    # Standard error log
#SBATCH --array=0-7                   # Array range

List=(1_A619 2_rel2 3_bif3 4_relk1 1_A619_2 2_rel2_2 3_bif3_2 4_relk1_2)
#List=(axillarybud1  axillarybud2  crownRoot1	crownRoot2)
#List=(Tassel_ba1 Tassel_bif1 Seedling_IAA2)

Path=/scratch/sb14489/3.scATAC_flo

#cd "$Path"/3.SortedBam
cd "$Path"/2.Mapped
#cd "$Path"/3.SortedBam_Class

module load SAMtools/1.10-iccifort-2019.5.281


samtools flagstat -@ 12 ./"${List[SLURM_ARRAY_TASK_ID]}"/outs/possorted_bam.bam > "${List[SLURM_ARRAY_TASK_ID]}"_MappedReads.txt
#samtools flagstat -@ 32 "${List[SLURM_ARRAY_TASK_ID]}"_Rmpcr.bam > "${List[SLURM_ARRAY_TASK_ID]}"_MappedReads_AfterPicard.txt


#samtools view -@ 32 -c -f 1 -F 12 "${List[SLURM_ARRAY_TASK_ID]}"_Sorted.bam > "${List[SLURM_ARRAY_TASK_ID]}"_MappedReads.txt
#samtools view -@ 32 -c -f 1 -F 12 "${List[SLURM_ARRAY_TASK_ID]}"_Rmpcr.bam > "${List[SLURM_ARRAY_TASK_ID]}"_MappedReads_AfterPicard.txt

#samtools view -@ 32 -c -F 260  "${List[SLURM_ARRAY_TASK_ID]}"_Sorted.bam > "${List[SLURM_ARRAY_TASK_ID]}"_MappedReads.txt
#samtools view -@ 32 -c -F 260  "${List[SLURM_ARRAY_TASK_ID]}"_Rmpcr.bam > "${List[SLURM_ARRAY_TASK_ID]}"_MappedReads_AfterPicard.txt
