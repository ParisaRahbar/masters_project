module load anaconda3/personal
source activate my_bowtie2_env

BAM=/rds/general/user/pr422/projects/epinott/live/user_analysed_data/Parisa/ADvas_PRR_20230504_PU1/HIS_bam


#FRIP score

## 3.5 featureCoutns

### covert BED (the peaks) to SAF
awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $4, $1, $2+1, $3, "."}' ${BAM}/../peaks/merged_PU1_HIS_KCZ0346_KCZ0347_KCZ0348_mus_macs2_001_peaks_peaks.narrowPeak > ${BAM}/../peaks/merged_PU1_HIS_KCZ0346_KCZ0347_KCZ0348_mus_macs2_001_peaks.narrowPeak.saf

### count
featureCounts -p -a ${BAM}/../peaks/merged_PU1_HIS_KCZ0346_KCZ0347_KCZ0348_mus_macs2_001_peaks.narrowPeak.saf -F SAF -o ${BAM}/../peaks/readCountInPeaks.txt ${BAM}/merged_PU1_HIS_KCZ0346_KCZ0347_KCZ0348_mus.bam

# 1. prepare
## convert BAM (BAM used to call peaks) to BED
bedtools bamtobed -i ${BAM}/merged_PU1_HIS_KCZ0346_KCZ0347_KCZ0348_mus.bam | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' > ${BAM}/merged_PU1_HIS_KCZ0346_KCZ0347_KCZ0348_mus.PE2SE.tagAlign

# 2. total reads in BAM/BED
samtools view -c ${BAM}/merged_PU1_HIS_KCZ0346_KCZ0347_KCZ0348_mus.bam

wc -l ${BAM}/merged_PU1_HIS_KCZ0346_KCZ0347_KCZ0348_mus.PE2SE.tagAlign

# 3. count reads in peak regions

## 3.1 tagAlign, intersectBed -a tagAlign -b bed
time bedtools sort -i ${BAM}/../peaks/merged_PU1_HIS_KCZ0346_KCZ0347_KCZ0348__mus_macs2_001_peaks.narrowPeak | bedtools merge -i stdin | bedtools intersect -u -a ${BAM}/merged_PU1_HIS_KCZ0346_KCZ0347_KCZ0348_mus.PE2SE.tagAlign -b stdin | wc -l