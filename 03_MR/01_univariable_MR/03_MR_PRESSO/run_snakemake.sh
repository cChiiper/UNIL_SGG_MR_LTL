#!/bin/bash

#SBATCH --job-name='MR-PRESSO'
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --partition urblauna
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=500MB

snakemake -j 30 --latency-wait=40 --cluster "sbatch -p urblauna --job-name='MR_PRESSO' --time=14:00:00 --nodes=1 --cpus-per-task=1 --mem 4GB"

