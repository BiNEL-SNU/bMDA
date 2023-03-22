#!/bin/bash

### Variable setup
groupID=190422
TUMOR_NAME=$groupID"_Tumor"
NORMAL_NAME=$groupID"_Normal"
GATK_NAME=201229-bMDA-BRCA-$groupID
GATK_SUFFIX=recal.snps.PASS.recode.vcf
VS_SUFFIX=varscan.snp.Somatic.hc.filtered  #최대한 많이 걸러낸 것임. 더 많은 variant를 원하면 loose하게 할 수 있음
MUTECT_SUFFIX=mutect1.final.vcf
LIST=list.file_$groupID #file containing list of file prefix
LIST_TUMOR=list.file_"$groupID"_Tumor #file containing list of file prefix (excluding matching normal pair)
EXOME_TARGET_BED=/home/share/reference/target/SS5.GRCh37.bed
ANNOVAR=/home/share/tools/annovar/ 
DATALOC=/home/jinhyun/NGS77/demux #location of bam raw data


### process GATK (input : $GATK_NAME.gatk.$GATK_SUFFIX file which is an output GATK UnifiedGenotyper)
vcf-subset $GATK_NAME.gatk.$GATK_SUFFIX -c $NORMAL_NAME -e > $NORMAL_NAME.$GATK_SUFFIX;
bedtools subtract -wa -u -header -a $GATK_NAME.gatk.$GATK_SUFFIX -b $NORMAL_NAME.$GATK_SUFFIX > $GATK_NAME.GATK.somatic.snp.vcf; 
cat $GATK_NAME.GATK.somatic.snp.vcf | grep -v ^# | wc -l #vcf header의 sample이름이 실제 파일 이름과 다를 경우 직접 바꿔준다.
while read line; do vcf-subset $GATK_NAME.GATK.somatic.snp.vcf -c $line -e > $GATK_NAME.Tumor.GATK.somatic.snp.$line.vcf; done < $LIST  #-e 옵션은 0/0 call이 아닌것만 남김

### Step 1. in-sample double calling
## inputs
## - GATK output from previous step
## - VARSCAN output ending with $VS_SUFFIX
## - MUTECT output ending with $MUTECT_SUFFIX

while read line; do 
    python ~/sourcecode/varscan.convert.vcf.py $line.$VS_SUFFIX > $line.$VS_SUFFIX.vcf
    cat $line.$VS_SUFFIX | awk '{OFS="\t"}{if ($1 != "chrom") print $1, ($2-1), $2}' > $line.$VS_SUFFIX.bed;

    bedtools intersect -a $GATK_NAME.Tumor.GATK.somatic.snp.$line.vcf -b $line.$VS_SUFFIX.vcf -header > $line.GATK.varscan.vcf;
    bedtools intersect -a $GATK_NAME.Tumor.GATK.somatic.snp.$line.vcf -b $line.$MUTECT_SUFFIX -header > $line.GATK.MuTect.vcf;
    bedtools intersect -a $line.$MUTECT_SUFFIX -b $line.$VS_SUFFIX.vcf -header > $line.MuTect.varscan.vcf;

    cat $line.GATK.varscan.vcf $line.GATK.MuTect.vcf $line.MuTect.varscan.vcf | grep -v ^# | awk '{OFS="\t"}{print $1, ($2-1),$2}' > $line.temp.bed;
    ~/sourcecode/sort.GRCh37.bed.sh $line.temp;
    bedtools merge -i $line.temp.bed > $line.double.called.bed; rm $line.temp.bed
done < $LIST_TUMOR

### Step 2. double called sites selection
INPUT=""; SF=double.called.bed; while read fname; do INPUT="$INPUT $fname.$SF"; done < $LIST_TUMOR
DCS=$GATK_NAME.double.called.sites.snp
cat $INPUT > temp.bed
~/sourcecode/sort.GRCh37.bed.sh temp
uniq temp.bed > $DCS.bed; ~/sourcecode/sort.GRCh37.bed.sh $DCS ; rm temp.bed; wc -l $DCS.bed; #total N sites 
bedtools intersect -a $DCS.bed -b $EXOME_TARGET_BED > $DCS.ontarget.bed; wc -l $DCS.ontarget.bed #N sites are on target

### Step 3. confident sites (inter sample N-called)
while read line; do cp $GATK_NAME.Tumor.GATK.somatic.snp.$line.vcf $line.GATK.vcf; done < $LIST_TUMOR
python ~/sourcecode/vcf.double.call.py $DCS.ontarget.bed $LIST_TUMOR confident.site 2
cat confident.site.variant  | grep -v , > confident.site.01.variant

### Step 4. extract raw data on confident sites (caustion : we need BAM raw data)
## example output (*.cs.snp.result) can be found in 'example'
cat $LIST | awk -v d=$DATALOC '{print d"/"$1"/"$1}' | grep -v $NORMAL_NAME > list.filepre
cat confident.site.01.variant | awk '{OFS="\t"}{print $1, $2, $2}' > confident.site.list
~/sourcecode/perform.bam.readcount.sh list.filepre confident.site.list 30 cs.snp 71 18 ##base quality > 30, mapping quality > 30
while read line; do python ~/sourcecode/perbase.parse.py confident.site.01 $line.cs.snp > $line.cs.snp.result; cp $line.cs.snp.result ./; done < list.filepre

### Step 5. Variant Call
## Either GATK, Varscan, MuTect, or Fisher's exect test (10^-3) PASSed -< Variant called
## (ref count + alt count < 5) -> NA 
## remaining: Variant not detected
python ~/sourcecode/vcf.call.by.rawdata.py $DCS.ontarget.bed $LIST_TUMOR confident.site 2
cat confident.site.filtered.variant |grep -v ^CHROM | awk '{OFS="\t"}{print $1, $2, $2, $3, $4}' > confident.site.filtered.avinput
INPUT=confident.site.filtered; 
$ANNOVAR/table_annovar.pl $INPUT.avinput $ANNOVAR/humandb/ -buildver hg19 -out $INPUT -remove -protocol refGene,1000g2015aug_all,avsnp150,cosmic70,icgc28,dbnsfp42c -operation g,f,f,f,f,f -nastring . 
