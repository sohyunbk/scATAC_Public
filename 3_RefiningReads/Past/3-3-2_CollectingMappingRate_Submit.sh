#!/bin/bash
#SBATCH --job-name=RunPythonForStat        # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=100gb                     # Job memory request
#SBATCH --time=3:00:00               # Time limit hrs:min:sec

python 3-3-2_CollectingMappingRate.py
