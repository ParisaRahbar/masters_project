#!/bin/bash

#
#PBS  -N Merging_peakcalling_frip_normalising
#
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=15:mem=60gb

module load anaconda3/personal
source activate my_bowtie2_env

export sample=('JVD0485_homo_P63_ERG_CUTTag_TF_WT_20230327_R1' 'JVD0488_homo_P94_ERG_CUTTag_TF_WT_20230310_R1' 'JVD0494_homo_UCSD_20221507_ERG_CUTTag_TF_WT_20230331_R1.')
BAM=/rds/general/user/pr422/projects/epinott/live/user_analysed_data/Parisa/human/ERG/TF/TF_bam

#Filtering out unmapped reads (-F 4) and low quality reads (-q 30)
for sample_type in "${sample[@]}"
do
samtools view -b -F 4 -q 30 ${BAM}/${sample_type}.target.linear_dedup.sorted.bam > ${BAM}/${sample_type}.target.linear_dedup.sorted.q30.bam 
done

#Mergin BAM files into one 
samtools merge ${BAM}/merged_ERG_TF_JVD0485_JVD0488_KCZ0455_homo.bam ${BAM}/*.q30.bam 

#Calling narrow peaks using macs2 
macs2 callpeak -t ${BAM}/merged_ERG_TF_JVD0485_JVD0488_KCZ0455_homo.bam \
 	-f BAM -g 2.9e+9 \
	-n merged_ERG_TF_JVD0485_JVD0488_KCZ0455_homo_macs2_001 \
	-q 0.01 \
	--outdir ${BAM}/../peaks 
#FRIP score

## 3.5 featureCoutns

### covert BED (the peaks) to SAF
awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $4, $1, $2+1, $3, "."}' ${BAM}/../peaks/merged_ERG_TF_JVD0485_JVD0488_KCZ0455_homo_macs2_001_peaks.narrowPeak > ${BAM}/../peaks/merged_ERG_TF_JVD0485_JVD0488_KCZ0455_homo_macs2_001_peaks.narrowPeak.saf


### count
featureCounts -p -a ${BAM}/../peaks/merged_ERG_TF_JVD0485_JVD0488_KCZ0455_homo_macs2_001_peaks.narrowPeak.saf -F SAF -o ${BAM}/../peaks/readCountInPeaks.txt ${BAM}/merged_ERG_TF_JVD0485_JVD0488_KCZ0455_homo.bam

# 1. prepare
## convert BAM (BAM used to call peaks) to BED
bedtools bamtobed -i ${BAM}/merged_ERG_TF_JVD0485_JVD0488_KCZ0455_homo.bam | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' > ${BAM}/merged_ERG_TF_JVD0485_JVD0488_KCZ0455_homo.PE2SE.tagAlign

# 2. total reads in BAM/BED
samtools view -c ${BAM}/merged_ERG_TF_JVD0485_JVD0488_KCZ0455_homo.bam

wc -l ${BAM}/merged_ERG_TF_JVD0485_JVD0488_KCZ0455_homo.PE2SE.tagAlign

# 3. count reads in peak regions

## 3.1 tagAlign, intersectBed -a tagAlign -b bed
time bedtools sort -i ${BAM}/../peaks/merged_ERG_TF_JVD0485_JVD0488_KCZ0455_homo_macs2_001_peaks.narrowPeak | bedtools merge -i stdin | bedtools intersect -u -a ${BAM}/merged_ERG_TF_JVD0485_JVD0488_KCZ0455_homo.PE2SE.tagAlign -b stdin | wc -l



#Making bigwig and normalising 

samtools index ${BAM}/merged_ERG_TF_JVD0485_JVD0488_KCZ0455_homo.bam

bamCoverage \
--bam ${BAM}/merged_ERG_TF_JVD0485_JVD0488_KCZ0455_homo.bam \
-o ${BAM}/merged_ERG_TF_JVD0485_JVD0488_KCZ0455_homo_linear_CPM.bam.bw \
--binSize 50 \
--normalizeUsing CPM \
--effectiveGenomeSize 2900000000 \
--minMappingQuality 30 


exit 0
