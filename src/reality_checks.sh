#!/bin/bash
#SBATCH -A zeller
#SBATCH -t 22:00:00
#SBATCH --mem 64G
#SBATCH -o /scratch/sim_bench/reality_checks/reality_%A.out
#SBATCH -e /scratch/sim_bench/reality_checks/reality_%A.err

module load Anaconda3
source activate /g/scb/zeller/jawirbel/software/envs/r_4.0.0/

cd ../
Rscript src/reality_checks_combination.R $1


