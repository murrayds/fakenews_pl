#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --job-name=Fake_news
#SBATCH --mem=100G
#SBATCH --partition=short


source activate fakenews
#echo $python 3.7.7

snakemake -j 1
