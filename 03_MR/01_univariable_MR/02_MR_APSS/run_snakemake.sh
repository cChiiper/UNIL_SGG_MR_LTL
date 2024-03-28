#!/bin/bash

#SBATCH --job-name='MR_APPS_launcher'
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --partition urblauna
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=500MB

snakemake -j 2 --latency-wait=40 --keep-going --cluster "sbatch -p urblauna --job-name='MR_APPS' --time={resources.time} --nodes={resources.nodes} --cpus-per-task={resources.cpus} --mem={resources.mem}"

