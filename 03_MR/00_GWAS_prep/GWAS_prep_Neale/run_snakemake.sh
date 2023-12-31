#!/bin/bash

#SBATCH --job-name='GWASprep'
#SBATCH --nodes=1
#SBATCH --time=5:00:00
#SBATCH --partition urblauna
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=2GB

snakemake -j 11 --cores=11 --cluster "sbatch -p urblauna --job-name='GWASprep' --time=1:00:00 --nodes=1 --cpus-per-task=4 --mem-per-cpu=15G"

