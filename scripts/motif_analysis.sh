#!/bin/bash

#
#PBS  -N motif_analysis
#
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=15:mem=60gb
#PBS -o Merging_peakcalling_frip_normalising.log

module load homer/4.11


findMotifsGenome.pl /rds/general/user/pr422/projects/epinott/live/user_analysed_data/Parisa/human/OLIG2/TF/peaks/KCZ0156_homo_P94_OLIG2_CUTTag_rIgG_50K_20221005_R1.macs2_summits.bed hg38 /rds/general/user/pr422/projects/epinott/live/user_analysed_data/Parisa/human/OLIG2/TF/motif_analysis -size 200 -mask -preparsedDir /rds/general/user/pr422/projects/epinott/live/user_analysed_data/Parisa/human/OLIG2/TF/motif_analysis


exit 0