#!bin/bash
set -e
set -u

date

SAMTOOLS="samtools"
REFERENCE=/home/tjamann/Documents/LRS/extdata/Zea_mays.AGPv3.27.dna.genome.fa
INTERVALS=/home/tjamann/Documents/LRS/extdata/IAD66395_209_Designed.interval_list
KNOWN=/home/tjamann/Documents/LRS/extdata/vcfsort.resorted.vcf
now=$(date +"%Y-%m-%d")
GATK="java -jar /home/tjamann/bin/GenomeAnalysisTK.jar"
PICARD="picard"
INTERVALS_BED=/home/tjamann/Documents/LRS/extdata/IAD66395_209_Designed.bed
DICT=/home/tjamann/Documents/LRS/extdata/Zea_mays.AGPv3.27.dna.genome.dict
HEADERFILE=/home/tjamann/Documents/LRS/extdata/combinedtargetmetrics.txt
KEEP=/home/tjamann/Documents/LRS/extdata/keep_seq.txt
HAPMAP=/home/tjamann/Documents/LRS/extdata/hapmap.sorted.vcf
outmode="EMIT_ALL_CONFIDENT_SITES"
emit_thresh=20	#Threshold for tagging possible variants
call_thresh=30	#Threshold for tagging _good_ variants
#unpack file and move to results
echo "unpack file"

for e in *.tar.bz2
do
  tar -xvjf $e --strip=1
done

LINES=20000
for f in *.fastq; do
  a=`cat "$f" | wc -l`;
  if [ "$a" -lt "$LINES" ]
  then
    rm -f "$f"
  fi
done

#number of reads that align to reference
#cat *.fastq > all.fastq
#bwa mem -t 4 /home/tjamann/Documents/LRS/extdata/Zea_mays.AGPv3.27.dna.genome.fa all.fastq > all.bwamem.sam
#samtools view -h -b -S all.bwamem.sam > all.bwamem.bam
#samtools sort all.bwamem.bam -o all.bwamem.sorted
#samtools index all.bwamem.sorted.bam
#samtools idxstats all.bwamem.sorted.bam > all.bwamem.bam.stats.txt
#bedtools coverage -a $INTERVALS_BED -b all.bwamem.sorted.bam > representationofamplicons.txt
#rm *.sam
#rm all.fastq
#rm all.bwamem.bam

#bwa mem alignment

echo "bwa mem alignment on all fastq files"

for f in *.fastq
do
  bwa mem -t 4 -R '@RG\tID:'$now'\tLB:P07\tSM:'$f'\tPL:IONTORRENT' $REFERENCE $f > $f.bwamem.sam
  samtools view -h -b -S $f.bwamem.sam > $f.bwamem.bam
  samtools view -b -F 4 $f.bwamem.bam > $f.bwamem.mapped.bam
  samtools sort $f.bwamem.mapped.bam -o $f.bwamem.mapped.sorted.bam
  samtools index $f.bwamem.mapped.sorted.bam
done

#echo #bowtie2 alignment
#for f in *.fastq
#do
#  bowtie2 -p 4 -x /media/sf_LRS/data/extdata/Zea_mays.AGPv3.27.dna.genome --sensitive-local $f > $f.bowtie.sam
#  samtools view -h -b -S $f.bowtie.sam > $f.bowtie.bam
#  samtools view -b -F 4 $f.bowtie.bam > $f.bowtie.mapped.bam
#  samtools sort $f.bowtie.mapped.bam $f.bowtie.mapped.sorted
#  samtools index $f.bowtie.mapped.sorted.bam
#done

echo "removing .mapped bam and .sam files"

rm *.mapped.bam
rm *.sam
rm *.fastq

#collecting sequence alignment metrics

echo "Assessing alignment" 

for BAM in *.mapped.sorted.bam; do picard CollectTargetedPcrMetrics I=$BAM O=${BAM%.bam}.targetmetrics.txt AMPLICON_INTERVALS=/home/tjamann/Documents/LRS/extdata/IAD66395_209_Designed.interval_list TARGET_INTERVALS=/home/tjamann/Documents/LRS/extdata/IAD66395_209_Designed.interval_list REFERENCE_SEQUENCE=/home/tjamann/Documents/LRS/extdata/Zea_mays.AGPv3.27.dna.genome.fa; done

for t in *.targetmetrics.txt
do 
  sed -e '1,7d' ${t} > ${t}.trimmed.txt
done

for t in *.trimmed.txt
do 
  sed -e '2,3d' ${t} > ${t}.trimmed.trimmed.txt
done

for i in *.trimmed.trimmed.txt
 do awk '{print FILENAME"\t"$0}' $i > $i.bk
done

rm *.trimmed.trimmed.txt
rm *.trimmed.txt

cat *.bk > out.txt

sed -e '1,6d' $HEADERFILE > header.txt
sed -e '2,4d' header.txt > newheader.txt
awk '{print "Filename""\t"$0}' newheader.txt > newnewheader.txt
rm newheader.txt
rm header.txt
cat newnewheader.txt out.txt > combinedtargetmetrics.txt
rm newnewheader.txt
rm *.bk
rm out.txt
rm *.sorted.targetmetrics.txt

find *.bwamem.mapped.sorted.bam -print > bwamem.all.bam.list
BAMLIST=bwamem.all.bam.list

$GATK \
-R $REFERENCE \
-T DepthOfCoverage \
-I $BAMLIST \
-L $INTERVALS \
-o DepthofCoverage

$GATK \
-R $REFERENCE \
-T FindCoveredIntervals \
-I $BAMLIST \
-L $INTERVALS \
-o FindCoveredIntervals

#picard CollectTargetedPcrMetrics I=all.bwamem.sorted.bam O=all.bwamem.sorted.bam.targetmetrics.txt AMPLICON_INTERVALS=$INTERVALS_BED TARGET_INTERVALS=$INTERVALS_BED REFERENCE_SEQUENCE=$REFERENCE
#rm all.bwamem.sorted.bam
#rm all.bwamem.sorted.bam.bai

## GATK Data Pre-Processing

# Step 1 - Local realignment around indels.
# Create a target list of intervals to be realigned.

echo "Creating a target list of intervals to be realigned...., local realignment, and indexing"

for B in *.bwamem.mapped.sorted.bam
do 
  $GATK \
  -T RealignerTargetCreator \
  -R $REFERENCE \
  -I $B \
  -o ${B%.bam}.target_intervals.list

  $GATK \
  -T IndelRealigner \
  -R $REFERENCE \
  -I $B \
  -targetIntervals "${B%.bam}.target_intervals.list" \
  -o ${B%.bam}.realigned_reads.bam

  $SAMTOOLS index "${B%.bam}.realigned_reads.bam"
done

find *.realigned_reads.bam -print > bwamem.realigned.bam.list
REALIGNEDBAMLIST=bwamem.realigned.bam.list

echo "calling variants"

outmode="EMIT_ALL_CONFIDENT_SITES"

for BAM in *.realigned_reads.bam
do
  $GATK \
  -T HaplotypeCaller \
  -R $REFERENCE \
  -I $BAM \
  -L $INTERVALS \
  -stand_call_conf 2.0 \
  -stand_emit_conf 1.0 \
  -out_mode $outmode \
  --emitRefConfidence GVCF \
  --variant_index_type LINEAR \
  --variant_index_parameter 128000 \
  -o $BAM.output.raw.snps.indels.g.vcf
done

echo "GATK Combine .g.vcf files"
# Run this after all of the variants have been called

find *.g.vcf -print > all.gvcf.list
GVCFLIST=all.gvcf.list

$GATK \
-T GenotypeGVCFs \
-R $REFERENCE \
--variant $GVCFLIST \
-stand_emit_conf $emit_thresh \
-stand_call_conf $call_thresh \
-o allsamples.raw.GATK.vcf \
-D $KNOWN
echo "outputting info about variants before filtering"



samtools mpileup -uf $REFERENCE --bam-list $REALIGNEDBAMLIST | bcftools call -vmO z -o allsamples.raw.mpileup.vcf.gz
tabix -p vcf allsamples.raw.mpileup.vcf.gz
###
#this script will filter variants after they are called using mpileup and GATK

#break vcf file so that multiallelic variants are in multiple lines instead of one single line. note- vcflib doesn't split it properly, so I use bcftools
#bgzip allsamples.raw.GATK.vcf
bgzip allsamples.raw.GATK.vcf
tabix -p vcf allsamples.raw.GATK.vcf.gz
bcftools norm -m-both -o GATK.break.step1.vcf allsamples.raw.GATK.vcf.gz
bcftools norm -f $REFERENCE -o allsamples.raw.GATK.break.vcf GATK.break.step1.vcf
rm GATK.break.step1.vcf

echo "filtering variants"
#filter here using GATK for quality scores
$GATK \
-R $REFERENCE \
-T VariantFiltration \
-o allsamplesafterfiltering.GATK.break.vcf \
--variant allsamples.raw.GATK.break.vcf \
--filterExpression "QD < 5.0" \
--filterName QDFilter \
--filterExpression "QUAL < 30.0" \
--filterName QUALFilter \
--genotypeFilterExpression "DP < 12" \
--genotypeFilterName "genotype_depth" \
--setFilteredGtToNocall
$GATK \
-R $REFERENCE \
-T VariantFiltration \
-o allsamples.mpileup.vcf \
--variant allsamples.raw.mpileup.vcf.gz \
--filterExpression "QUAL < 30.0" \
--filterName QUALFilter \
--genotypeFilterExpression "DP < 12" \
--genotypeFilterName "genotype_depth" \
--setFilteredGtToNocall
rm allsamples.raw.GATK.break.vcf

#use vcftools to filter for minor allele frequency and to filter out indels so only SNPs remain. also filter the hapmap file the same (except for maf) for later use
vcftools --vcf allsamplesafterfiltering.GATK.break.vcf --remove-filtered-all --remove-indels --maf 0.02 --recode --out final.GATK.break.vcf
vcftools --vcf allsamples.mpileup.vcf --remove-filtered-all --remove-indels --recode --maf 0.02 --out final.mpileup
vcftools --vcf $HAPMAP --remove-filtered-all --remove-indels --recode --out final.hapmap
rm allsamplesafterfiltering.GATK.break.vcf
rm allsamples.mpileup.vcf

#filter mpileup so that only targeted regions are kept. GATK is already filtered this way (only regions of interest are used to call variants)
bedtools intersect -a final.mpileup.recode.vcf -b $INTERVALS_BED -header > allsamplesafterfiltering.mpileup.vcf

echo "outputting info about variants after filtering"

#get information about variants that are left.
$GATK \
-R $REFERENCE \
-T VariantEval \
-L $INTERVALS \
-o mpileup.eval.gatkreport \
--eval: allsamplesafterfiltering.mpileup.vcf \
-D $KNOWN

$GATK \
-R $REFERENCE \
-T VariantEval \
-L $INTERVALS \
-o GATK.eval.gatkreport \
--eval: final.GATK.break.vcf.recode.vcf \
-D $KNOWN

#compress so that can compare vcf files from different programs using vcf-compare. Used vcftools because GATK takes longer and is having trouble with the GATK file (reference sorting)
bgzip allsamplesafterfiltering.mpileup.vcf
tabix -p vcf allsamplesafterfiltering.mpileup.vcf.gz
bgzip final.GATK.break.vcf.recode.vcf
tabix -p vcf final.GATK.break.vcf.recode.vcf.gz
bgzip final.hapmap.recode.vcf
tabix -p vcf final.hapmap.recode.vcf.gz
vcf-compare -a allsamplesafterfiltering.mpileup.vcf.gz final.GATK.break.vcf.recode.vcf.gz final.hapmap.recode.vcf.gz > compare.txt

echo "combining files with positive controls"
#filter so that only inbreds and checks remain and then combine into a vcf file.
vcftools --gzvcf allsamplesafterfiltering.mpileup.vcf.gz --keep $KEEP --remove-filtered-all --recode --out mpileup.filteredindiv
vcftools --gzvcf final.GATK.break.vcf.recode.vcf.gz --keep $KEEP --remove-filtered-all --recode --out GATK.filteredindiv
bgzip mpileup.filteredindiv.recode.vcf
bgzip GATK.filteredindiv.recode.vcf
tabix -p vcf mpileup.filteredindiv.recode.vcf.gz
tabix -p vcf GATK.filteredindiv.recode.vcf.gz
vcf-merge mpileup.filteredindiv.recode.vcf.gz GATK.filteredindiv.recode.vcf.gz final.hapmap.recode.vcf.gz > mpile.GATK.true.vcf
vcftools --vcf mpile.GATK.true.vcf --012

rm *.idx

#analyze discovered variants
snpEff -v -csvStats snpeffbed.txt AGPv3.27 -i bed $INTERVALS_BED -s summaryintervals > bed_targets.vcf
snpEff -v -csvStats snpeffsamtools.txt AGPv3.27 allsamplesafterfiltering.mpileup.vcf -s summarysamtools > samtools_annotation.vcf
snpEff -v -csvStats snpeffGATK.txt AGPv3.27 final.GATK.break.vcf.recode.vcf -s summaryGATK > GATK_annotation.vcf
snpEff -v -csvStats snpeffhapmap.txt AGPv3.27 final.hapmap.recode.vcf -s summaryhapmap > hapmap_annotation.vcf


#for input to tassel GATK snps for tree construction

bgzip -d allsamples.raw.GATK.vcf.gz

$GATK \
-R $REFERENCE \
-T VariantFiltration \
-o allsamplesafterfiltering.GATK.vcf \
--variant allsamples.raw.GATK.vcf \
--filterExpression "QD < 5.0" \
--filterName QDFilter \
--filterExpression "QUAL < 30.0" \
--filterName QUALFilter \
--genotypeFilterExpression "DP < 12" \
--genotypeFilterName "genotype_depth" \
--setFilteredGtToNocall

vcftools --vcf allsamplesafterfiltering.GATK.vcf --remove-indels --remove-filtered-all --maf 0.02 --recode --recode-INFO-all --out SNPs_only
bcftools annotate -x INFO,^FORMAT/GT SNPs_only.recode.vcf > trythis.vcf
sed '/*/d' trythis.vcf > trythissed.vcf
/home/tjamann/bin/tassel-5-standalone/run_pipeline.pl -SortGenotypeFilePlugin -inputFile trythissed.vcf -outputFile sortedtassel.vcf -fileType VCF
date