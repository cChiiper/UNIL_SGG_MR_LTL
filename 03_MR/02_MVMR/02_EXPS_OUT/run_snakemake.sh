#!/bin/bash

#SBATCH --job-name='SN-MVMR'
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --partition urblauna
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=500MB

snakemake -j 3 --latency-wait=40 --cluster "sbatch -p urblauna --job-name='MVMR_test' --time=00:15:00 --nodes=1 --cpus-per-task=1 --mem 20GB"

