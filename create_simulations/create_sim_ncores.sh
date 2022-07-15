#!/bin/bash
#SBATCH -A zeller
#SBATCH -t 3-00:00:00
#SBATCH --mem 64G
#SBATCH --cpus-per-task=64
#SBATCH -o /scratch/jawirbel/simulation_benchmark/github/create_sim_%A.out
#SBATCH -e /scratch/jawirbel/simulation_benchmark/github/create_sim_%A.err

module load R/4.0.0-foss-2020a

cd ../

Rscript ./create_simulations/$1
