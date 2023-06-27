#!/bin/bash

#
#PBS  -N Merging_peakcalling_frip_normalising
#
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=15:mem=60gb

module load anaconda3/personal
source activate my_bowtie2_env

export sample=('KCZ0660_homo_UCSD_20220803_PU1_nanoCUTTag_rb_H3K27ac_WT_20230523_R1')
BAM=/rds/general/user/pr422/projects/epinott/live/user_analysed_data/Parisa/human_nano/PU1/bam_files

#Filtering out unmapped reads (-F 4) and low quality reads (-q 30)
for sample_type in "${sample[@]}"
do
samtools view -b -F 4 -q 30 ${BAM}/${sample_type}.target.linear_dedup.sorted.bam > ${BAM}/${sample_type}.target.linear_dedup.sorted.q30.bam 
done