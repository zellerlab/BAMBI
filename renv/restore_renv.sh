#!/bin/bash
#SBATCH -A zeller
#SBATCH -t 3-00:00:00
#SBATCH --mem 64G
#SBATCH -o /scratch/sim_bench/renv.out
#SBATCH -e /scratch/sim_bench/renv.err

module load R/4.0.0-foss-2020a

cd ../
Rscript ./renv/restore_renv.R
