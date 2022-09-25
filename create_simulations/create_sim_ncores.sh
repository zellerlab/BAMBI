#!/bin/bash
#SBATCH -A zeller
#SBATCH -t 3-00:00:00
#SBATCH --mem 64G
#SBATCH --cpus-per-task=64
#SBATCH -o /scratch/jawirbel/simulation_benchmark/github/create_sim_%A.out
#SBATCH -e /scratch/jawirbel/simulation_benchmark/github/create_sim_%A.err

module load Anaconda3
source activate /g/scb/zeller/jawirbel/software/envs/r_4.0.0

cd ../

Rscript ./create_simulations/$1
