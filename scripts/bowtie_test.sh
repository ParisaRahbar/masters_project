#!/bin/bash

#
#PBS  -N BOWTIETEST
#
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=10:mem=15gb

module load anaconda3/personal
source activate my_bowtie2_env

projPath="/rds/general/user/pr422/home/projPath"
ref="/rds/general/user/pr422/home/projPath/alignment/mouse/mm10_index/mm10"
read1="/rds/general/user/pr422/home/projPath/trimmed/lhx2/IGF123984_val_1.fq.gz"
read2="/rds/general/user/pr422/home/projPath/trimmed/lhx2/IGF123984_val_2.fq.gz"
export TMPDIR="/rds/general/user/pr422/home/projPath/Cut_Tag_analysis/tmp_dir"
cores=8
bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${ref} -1 ${read1} -2 ${read2} -S ${projPath}/alignment/sam/IGF123984_bowtie2.sam &> ${projPath}/alignment/sam/bowtie2_summary/lhx2/IGF123984_bowtie2.txt

exit 0