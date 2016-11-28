#compress so that can compare vcf files from different programs using vcf-compare. Used vcftools because GATK takes longer and is having trouble with the GATK file (reference sorting)
bgzip allsamplesafterfiltering.mpileup.vcf
tabix -p vcf allsamplesafterfiltering.mpileup.vcf.gz
bgzip final.GATK.break.vcf.recode.vcf
tabix -p vcf final.GATK.break.vcf.recode.vcf.gz
bgzip final.hapmap.recode.vcf
tabix -p vcf final.hapmap.recode.vcf.gz
vcf-compare -a final.GATK.break.vcf.recode.vcf.gz final.hapmap.recode.vcf.gz > compare.txt

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



date