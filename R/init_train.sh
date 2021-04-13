#!/bin/bash
#SBATCH --ntasks-per-node=32
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=100M
#SBATCH --time=0:01:00
#SBATCH --job-name=job2
#SBATCH --error=job2.out
#SBATCH --output=job2.out

source /etc/bashrc
module load R/4.0.0
cd /cluster/home/$USER/text_features/R
Rscript init_train.R CCLE rf regression