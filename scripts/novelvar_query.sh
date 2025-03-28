#!/bin/sh

#$1 and $2 should be passed when calling the script. 
#$1 is the vcf file while $2 is the output file

#Extract novel variants from the data, those that have no rsid in the dbSNPs (annotated by '.' in chr position).
#bcftools view -n $1 -Oz -o $2novel_variants.vcf.gz
#tabix -p vcf $2novel_variants.vcf.gz

bcftools view -S ^samplestoremove.list -i 'MAF>0.01' $2novel_variants.vcf.gz -Oz -o $2novel_variants_common.vcf.gz
tabix -p vcf $2novel_variants_common.vcf.gz

#Annotated the novel variants common (MAF>0.01)
java -Xmx20g -jar ~/Softwares/snpEff/snpEff.jar -v hg38 data/novel_variants_common.vcf.gz > output/novel_variants_common.annot.vcf
bgzip -c output/novel_variants_common.annot.vcf > output/novel_variants_common.annot.vcf.gz

#Produce summary stats of the output/novel_variants_common.annot.vcf using bcftools and plotted it
bcftools stats -s - output/novel_variants_common.annot.vcf > output/novel_variants_common.annot.stats
plot-vcfstats -p plot/ output/novel_variants_common.annot.stats

#Filtered the common novel variants with AF>=0.9 to identify alleles that are highly frequent in the Filipino population
bcftools view -i 'AF>=0.99' output/novel_variants_common.annot.vcf -Oz -o output/novel_variants_common_af99.vcf.gz
bcftools view -i 'AF>=0.98' output/novel_variants_common.annot.vcf.gz -Oz -o output/novel_variants_common_af98.vcf.gz

#Run the ROH
bcftools roh --AF-tag AF ../output/novel_variants_common.annot.vcf -o ../output/novel_variants_common.annot.roh.txt

Stats	Var_Count
output/novel_variants_common_af90.vcf.gz 415495	

#Annotate using SIFT
java -jar ~/Softwares/SIFT4G_Annotator.jar -c -i ~/RawData/Fgrp02/GRCh38/20240124_genotypes.FINAL.vcf -d resources/GRCh38.83.chr/ -t -r output/sift/rawdata.sift
java -jar ~/Softwares/SIFT4G_Annotator.jar -c -i output/novel_variants_common.annot.vcf -d resources/GRCh38.83.chr/ -t -r output/sift/novel_variants_common.sift

# Annotate using polyphen
# Downloaded polyphen resource whess
wget http://genetics.bwh.harvard.edu/downloads/pph2/whess/polyphen-2.2.2-whess-2011_12.sqlite.bz2
bunzip2 polyphen-2.2.2-whess-2011_12.sqlite.bz2
