#!/bin/bash
#SBATCH -A zeller
#SBATCH -t 3-00:00:00
#SBATCH --mem 64G
#SBATCH -o /scratch/sim_bench/run_tests/test_%A.out
#SBATCH -e /scratch/sim_bench/run_tests/test_%A.err

module load R/4.0.0-foss-2020a

cd ../
Rscript ./src/run_tests.R $1 $2
