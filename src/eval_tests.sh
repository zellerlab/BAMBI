#!/bin/bash
#SBATCH -A zeller
#SBATCH -t 3-00:00:00
#SBATCH --mem 8G
#SBATCH -n 1
#SBATCH -o /scratch/jawirbel/sim_bench/eval/eval_%A.out
#SBATCH -e /scratch/jawirbel/sim_bench/eval/eval_%A.err

module load R/4.0.0-foss-2020a

cd ../
Rscript src/eval_tests.R $1 $2


