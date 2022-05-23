#!/bin/bash
#SBATCH -A zeller
#SBATCH -t 22:00:00
#SBATCH --mem 64G
#SBATCH -o /scratch/sim_bench/reality_checks/reality_%A.out
#SBATCH -e /scratch/sim_bench/reality_checks/reality_%A.err

module load R/4.0.0-foss-2020a

cd ../
Rscript src/reality_checks_combination.R $1


