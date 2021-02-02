#! /bin/bash -l
#$ -j y
#$ -m ea
#$ -l buyin
#$ -l h_rt=12:00:00
#$ -pe omp 4
#$ -l mem_per_core=16G
#$ -N nullNimbleTest
#$ -o ../Output_Files/RCNfx/

module load R/4.0.2
Rscript 03_RandomWalk.R
