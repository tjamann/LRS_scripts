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

echo "bwa mem alignment on all fastq files"

for f in *.fastq
do
  bwa mem -t 4 -R '@RG\tID:'$now'\tLB:P07\tSM:'$f'\tPL:IONTORRENT' $REFERENCE $f > $f.bwamem.sam
  samtools view -h -b -S $f.bwamem.sam > $f.bwamem.bam
  samtools view -b -F 4 $f.bwamem.bam > $f.bwamem.mapped.bam
  samtools sort $f.bwamem.mapped.bam -o $f.bwamem.mapped.sorted.bam
  samtools index $f.bwamem.mapped.sorted.bam
done

echo "removing .mapped bam and .sam files"

rm *.mapped.bam
rm *.sam
rm *.fastq

find *.bwamem.mapped.sorted.bam -print > bwamem.all.bam.list
BAMLIST=bwamem.all.bam.list

$GATK \
-R $REFERENCE \
-T DepthOfCoverage \
-I $BAMLIST \
-L $INTERVALS \
-o DepthofCoverage

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

###
#this script will filter variants after they are called using GATK

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

rm allsamples.raw.GATK.break.vcf

#use vcftools to filter for minor allele frequency and to filter out indels so only SNPs remain. also filter the hapmap file the same (except for maf) for later use
vcftools --vcf allsamplesafterfiltering.GATK.break.vcf --remove-filtered-all --remove-indels --maf 0.001 --recode --out final.GATK.break.vcf

rm allsamplesafterfiltering.GATK.break.vcf

echo "preparing file so tassel can read it"

#for file so that tassel can read it

$GATK \
-R $REFERENCE \
-T VariantFiltration \
-o final.filtered.GATK.vcf \
--variant final.GATK.break.vcf.recode.vcf \
--filterExpression "QD < 5.0" \
--filterName QDFilter \
--filterExpression "QUAL < 30.0" \
--filterName QUALFilter \
--genotypeFilterExpression "DP < 12" \
--genotypeFilterName "genotype_depth" \
--setFilteredGtToNocall

vcftools --vcf final.filtered.GATK.vcf --remove-indels --remove-filtered-all --maf 0.02 --recode --recode-INFO-all --out SNPs_only
bcftools annotate -x INFO,^FORMAT/GT SNPs_only.recode.vcf > trythis.vcf
sed '/*/d' trythis.vcf > trythissed.vcf
/home/tjamann/bin/tassel-5-standalone/run_pipeline.pl -SortGenotypeFilePlugin -inputFile trythissed.vcf -outputFile sortedtassel.vcf -fileType VCF
rm trythis.vcf

