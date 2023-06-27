#!/bin/bash

#
#PBS  -N Merging_peakcalling_frip_normalising
#
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=20:mem=30gb

module load anaconda3/personal
source activate base 

export sample=('KCZ0353_mus_CR017_ERGhi_CUTTag_rIgG_linear_20221122_R1' 'KCZ0444_mus_CR035_ERG_CUTTag_TF_WT_20230306_R1' 'KCZ0453_mus_CR038_ERG_CUTTag_TF_WT_20230306_R1')
BAM=/rds/general/user/pr422/projects/epinott/live/user_analysed_data/Parisa/ADvas_PRR_20230504_ERG/TF_bam

module load samtools

#Filtering out unmapped reads (-F 4) and low quality reads (-q 30)
for sample_type in "${sample[@]}"
do
samtools view -b -F 4 -q 30 ${BAM}/${sample_type}.target.dedup.sorted.bam > ${BAM}/${sample_type}.target.dedup.sorted.q30.bam 
done

#Mergin BAM files into one 
samtools merge ${BAM}/merged_ERG_TF_KCZ0353_KCZ0347_KCZ0453_mus.bam ${BAM}/*.q30.bam 

source activate my_bowtie2_env

#Calling narrow peaks using macs2 
macs2 callpeak -t ${BAM}/merged_ERG_TF_KCZ0353_KCZ0347_KCZ0453_mus.bam \
 	-f BAM -g 1.87e+9 \
	-n merged_ERG_TF_KCZ0353_KCZ0347_KCZ0453_mus_macs2_001_peaks \
	-q 0.01 \
	--outdir ${BAM}/../peaks 

#FRIP score

## 3.5 featureCoutns
module load subread

### covert BED (the peaks) to SAF
awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $4, $1, $2+1, $3, "."}' ${BAM}/../peaks/merged_ERG_TF_KCZ0353_KCZ0347_KCZ0453_mus_macs2_001_peaks.narrowPeak > ${BAM}/../peaks/merged_ERG_TF_KCZ0353_KCZ0347_KCZ0453_mus_macs2_001_peaks.narrowPeak.saf

module load bedtools

### count
featureCounts -p -a ${BAM}/../peaks/merged_ERG_TF_KCZ0353_KCZ0347_KCZ0453_mus_macs2_001_peaks.narrowPeak.saf -F SAF -o ${BAM}/../peaks/readCountInPeaks.txt ${BAM}/merged_ERG_TF_KCZ0353_KCZ0347_KCZ0453_mus.bam

# 1. prepare
## convert BAM (BAM used to call peaks) to BED
bedtools bamtobed -i ${BAM}/merged_ERG_TF_KCZ0353_KCZ0347_KCZ0453_mus.bam | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' > ${BAM}/merged_ERG_TF_KCZ0353_KCZ0347_KCZ0453_mus.PE2SE.tagAlign

# 2. total reads in BAM/BED
samtools view -c ${BAM}/merged_ERG_TF_KCZ0353_KCZ0347_KCZ0453_mus.bam

wc -l ${BAM}/merged_ERG_TF_KCZ0353_KCZ0347_KCZ0453_mus.PE2SE.tagAlign

# 3. count reads in peak regions

## 3.1 tagAlign, intersectBed -a tagAlign -b bed
time bedtools sort -i ${BAM}/../peaks/merged_ERG_TF_KCZ0353_KCZ0347_KCZ0453_mus_macs2_001_peaks.narrowPeak | bedtools merge -i stdin | bedtools intersect -u -a ${BAM}/merged_ERG_TF_KCZ0353_KCZ0347_KCZ0453_mus.PE2SE.tagAlign -b stdin | wc -l



#Making bigwig and normalising 

samtools index ${BAM}/merged_ERG_TF_KCZ0353_KCZ0347_KCZ0453_mus.bam

bamCoverage \
--bam ${BAM}/merged_ERG_TF_KCZ0353_KCZ0347_KCZ0453_mus.bam \
-o ${BAM}/merged_ERG_TF_KCZ0353_KCZ0347_KCZ0453_mus_CPM.bw \
--binSize 50 \
--normalizeUsing CPM \
--effectiveGenomeSize 1870000000 \
--minMappingQuality 30 


exit 0

