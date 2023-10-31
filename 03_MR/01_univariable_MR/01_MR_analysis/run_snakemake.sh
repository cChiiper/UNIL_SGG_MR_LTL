#!/bin/bash

#SBATCH --job-name='TwoSampleMR'
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --partition urblauna
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=2GB

snakemake -j 2 --latency-wait=40 --cluster "sbatch -p urblauna --job-name='TwoSampleMR' --time=1:00:00 --nodes=1 --cpus-per-task=1 --mem-per-cpu=10G"

