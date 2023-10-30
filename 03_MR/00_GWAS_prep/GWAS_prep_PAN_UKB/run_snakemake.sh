#!/bin/bash

#SBATCH --job-name='sn_GWAS'
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --partition normal
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=2GB

snakemake -j 2 --cluster "sbatch -p normal --job-name='GWASprep' --time=2:00:00 --nodes=1 --cpus-per-task=1 --mem-per-cpu=40G"

