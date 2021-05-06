#! /bin/bash -l
#$ -j y
#$ -m ea
#$ -l buyin
#$ -l h_rt=72:00:00
#$ -pe omp 8
#$ -l mem_per_core=18G
#$ -o ../Output_Files/RCNfx/

module load R/4.0.2

# #$ -N plotEffects
# Rscript 04.1_plotEffects.R

# #$ -N yearEffects
# Rscript 04.2_yearEffects.R

# #$ -N cumGDDDiff
# Rscript 04.3_indModelSetUp.R

# #$ -N smam
# Rscript 04.4_indSmamPrevYear.R

#$ -N obsPoisTempSurv
Rscript 04.5.1_groupEffectsPois.R
