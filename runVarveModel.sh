#!/bin/bash
 
#SBATCH --job-name="Varve Model"
#SBATCH --output="/home/npm4/VarveModel/vm.out"
#SBATCH --workdir="/home/npm4/VarveModel/"
#SBATCH -n 1
#SBATCH --time=6000:00
#SBATCH --cpus-per-task=1 
#SBATCH --mail-type=ALL
#SBATCH --mem=64000
module load R

srun Rscript varveModel.R
