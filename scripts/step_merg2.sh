#!/bin/bash

#
#PBS  -N Merging_peakcalling_frip_normalising
#
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=15:mem=60gb

module load anaconda3/personal
source activate my_bowtie2_env

export sample=('KCZ0641_homo_P94_ERGhi_CUTTag_H3K27ac_WT_20230509_R1')
BAM=/rds/general/user/pr422/projects/epinott/live/user_analysed_data/Parisa/human_cut/ERG/HIS/bam_files_peaks/ERGhigh


#Calling narrow peaks using macs2 
macs2 callpeak -t ${BAM}/${sample_type}.target.linear_dedup.sorted.q30.bam  \
 	-f BAM -g 2.9e+9 \
	-n merged_ERGhigh_HIS_KCZ0641_homo_macs2_001 \
	-q 0.01 \
	--outdir ${BAM}/../peaks 


