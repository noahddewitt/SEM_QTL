#!/bin/bash

#SBATCH -p medium 
#SBATCH --job-name="Kin19_graphIter"
#SBATCH -N 1 
#SBATCH -n 72 
#SBATCH -t 168:00:00
#SBATCH --mail-user=nddewitt@ncsu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o "stdout.%j.%N"

module load r

Rscript graphConstruction.R baseGraph 72

