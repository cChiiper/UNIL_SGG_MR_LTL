#!/bin/bash

#SBATCH --job-name='SENS-Steiger'
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --partition urblauna
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=500MB

snakemake -j 50 --latency-wait=40 --cluster "sbatch -p urblauna --job-name='SENS-Steiger' --time=02:00:00 --nodes=1 --cpus-per-task=1 --mem 16GB"

