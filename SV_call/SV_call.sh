#!/bin/bash

### variable setup
groupID=190422
DATALOC=/home/jinhyun/NGS73/demux/CNA/SVanalysis
SV_VCF_LOC=/home/kyungsub/bMDA
LIST_TUMOR=../list.file_"$groupID"_Tumor
TUMOR_BULK_BAM=/home/jinhyun/NGS68/WGS/resorted/190422_Tumor.resorted.bam

### setup for filtering out low mapping quality regions
samtools view -H $TUMOR_BULK_BAM > lowMQ.sam
samtools view $TUMOR_BULK_BAM | awk '$5<20 {print $0}' >>  lowMQ.sam
samtools view -S -b -h lowMQ.sam > lowMQ.bam
samtools depth lowMQ.bam >  lowMQ.cov
SURVIVOR bincov lowMQ.cov 10 2 > lowMQ.bed

### high confident somatic SVs across all samples
mkdir delly; mkdir gridss; mkdir manta
rm input_files
while read line; do #merge by caller for each samples
echo processing $line;
zcat $SV_VCF_LOC/gridss/"$line"_high_confidence_somatic.vcf.bgz | awk '{OFS = "\t"}{if ($0 !~ /^##/) {t = $10; $10 = $11; $11 = t}; print;}' > gridss/$line.vcf #"$line"_high_and_low_confidence_somatic.vcf.bgz
zcat $SV_VCF_LOC/manta/"$line"/results/variants/somaticSV.vcf.gz | awk '{OFS = "\t"}{if ($0 !~ /^##/) {t = $10; $10 = $11; $11 = t}; print;}' > manta/$line.vcf #| awk '{if (($0 ~ /^#/) || ($8 !~ /IMPRECISE/)) print $0;}' |
echo manta/$line.vcf > $line.input_files; 
echo gridss/$line.vcf >> $line.input_files;
SURVIVOR merge $line.input_files 1000 1 1 1 0 30 "$line"_merged.vcf
echo "$line"_merged.vcf >> input_files
done < $LIST_TUMOR

while read line; do echo processing $line; #check headers (only 1st samples are merged by SURVIVOR)
cat manta/$line.vcf | grep ^# | tail -n 1; cat delly/$line.vcf | grep ^# | tail -n 1; cat gridss/$line.vcf  | grep ^# | tail -n 1; echo ""
done < $LIST_TUMOR

SURVIVOR merge input_files 1000 2 1 1 0 30 190422_SV_MANTAandGRIDSS_merged.vcf

INPUT=190422_SV_MANTAandGRIDSS_merged
SURVIVOR filter $INPUT.vcf lowMQ.bed 50 -1 0.01 -1 $INPUT.filtered.vcf #RE flag is not found in VCF

SURVIVOR stats $INPUT.vcf -1 -1 -1 SV_summary
SURVIVOR stats $INPUT.filtered.vcf -1 -1 -1 SV_summary

python ~/sourcecode/SVvcf_analysis.py $INPUT.filtered.vcf $INPUT.filtered.variant


### high confident somatic SVs supported by bulk
rm input_files
while read line; do #merge by caller for each samples
echo processing $line;
cat $SV_VCF_LOC/delly/$line.vcf > delly/$line.vcf
zcat $SV_VCF_LOC/gridss/"$line"_high_and_low_confidence_somatic.vcf.bgz | awk '{OFS = "\t"}{if ($0 !~ /^##/) {t = $10; $10 = $11; $11 = t}; print;}' > gridss/$line.vcf #"$line"_high_and_low_confidence_somatic.vcf.bgz
zcat $SV_VCF_LOC/manta/"$line"/results/variants/somaticSV.vcf.gz | awk '{OFS = "\t"}{if ($0 !~ /^##/) {t = $10; $10 = $11; $11 = t}; print;}' | awk '{if (($0 ~ /^#/) || ($8 !~ /IMPRECISE/)) print $0;}' > manta/$line.vcf
echo manta/$line.vcf > $line.input_files; 
echo gridss/$line.vcf >> $line.input_files;
echo delly/$line.vcf >> $line.input_files;  
SURVIVOR merge $line.input_files 1000 1 1 1 0 30 "$line"_merged.vcf
echo "$line"_merged.vcf >> input_files
done < $LIST_TUMOR

while read line; do echo processing $line; #check headers (only 1st samples are merged by SURVIVOR)
cat manta/$line.vcf | grep ^# | tail -n 1; cat delly/$line.vcf | grep ^# | tail -n 1; cat gridss/$line.vcf  | grep ^# | tail -n 1; echo ""
done < $LIST_TUMOR

SURVIVOR merge input_files 1000 1 1 1 0 30 "$groupID"_SV_all_merged.vcf

cat "$groupID"_SV_all_merged.vcf | awk '{split($(NF), gt, ":"); split($8, info, ";"); split(info[1], SUPP, "="); if (($0 ~ /^#/) || ((gt[2] != "NaN") && (SUPP[2] >= 2))) print $0}' > "$groupID"_SV_all_merged_lowConf_bulkFilter.vcf


INPUT="$groupID"_SV_all_merged_lowConf_bulkFilter
SURVIVOR filter $INPUT.vcf lowMQ.bed 50 -1 0.01 -1 $INPUT.filtered.vcf #RE flag is not found in VCF

SURVIVOR stats $INPUT.vcf -1 -1 -1 SV_summary
SURVIVOR stats $INPUT.filtered.vcf -1 -1 -1 SV_summary

python ~/sourcecode/SVvcf_analysis.py $INPUT.filtered.vcf $INPUT.filtered.variant