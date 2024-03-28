#!/bin/bash

#SBATCH --job-name='MR_Pipeline'
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --partition urblauna
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=500MB

snakemake -j 151 --latency-wait=40 --keep-going --cluster "sbatch -p urblauna --job-name='MR_LTL' --time={resources.time} --nodes={resources.nodes} --cpus-per-task={resources.cpus} --mem={resources.mem}"

