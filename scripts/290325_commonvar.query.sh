#!/bin/sh

# Ensure correct argument usage
if [ "$1" != "--raw" ] || [ "$3" != "--out" ]; then
	    echo "Usage: $0 --raw <input.vcf.gz> --out <OutputDirectory>"
	        exit 1
fi

# Command-line input arguments
raw="$2"	#Path to the raw vcf (compressed) dataset
out="$4"	#Path to the directory for the output file

echo "Investigating All Variants"

echo "===Step 1: Removing related individuals, low quality variants, and samples with high-missing data==="
bcftools view -S ^samplestoremove.list -i 'MAF>0.01 && QUAL>30 && FMT/GQ>10 && FILTER="PASS"' ${raw} -Oz -o ${out}all_variants.common.vcf.gz
tabix -p vcf ${out}all_variants.common.vcf.gz
bgzip -d ${out}all_variants.common.vcf.gz > ${out}all_variants.common.vcf

echo "===Step 2: Annotating using SnpEff==="
#Annotate the common all variants (MAF>0.01) using SnpEff
java -Xmx20g -jar ~/Softwares/snpEff/snpEff.jar -v hg38 ${out}all_variants.common.vcf.gz -stats ${out}all_variants.common.snpEff > ${out}all_variants.common.snpEff.vcf
bgzip -c ${out}all_variants.common.snpEff.vcf > ${out}all_variants.common.snpEff.vcf.gz
rm ${out}all_variants.common.snpEff.vcf

echo "===Step 3: Annotating using SIFT==="
#Annotate using SIFT, single transcipt only
mkdir ${out}all_variants.common.sift
java -jar ~/Softwares/SIFT4G_Annotator.jar -c -i ${out}all_variants.common.vcf -d ../resources/GRCh38.83.chr/ -r ${out}all_variants.common.sift

sed -n '1p;/DELETERIOUS/p' ${out}all_variants.common.sift/all_variants.common_SIFTannotations.xls | grep -v "*WARNING!" > ${out}all_variants.common.sift_deleterious.txt
sed -n '1p;/all/p' ${out}all_variants.common.sift_deleterious.txt > ${out}all_variants.common.sift_deleterious_all.txt
sed -n '1p;/all/!p' ${out}all_variants.common.sift_deleterious.txt > ${out}all_variants.common.sift_deleterious_dbSNP.txt

count_del=$(awk -F '\t' '{print $6}' ${out}all_variants.common.sift_deleterious.txt | sort -u | wc -l)
count_del_all=$(awk -F '\t' '{print $6}' ${out}all_variants.common.sift_deleterious_all.txt | sort -u | wc -l)
count_del_dbSNP=$(awk -F '\t' '{print $6}' ${out}all_variants.common.sift_deleterious_dbSNP.txt | sort -u | wc -l)

#need to subtract 1 to remove the header line count since we counted it including the header line
echo "The input callset has $(expr $count_del - 1) deleterious variant, of which $(expr $count_del_all - 1) are found to be all to the dataset and $(expr $count_del_dbSNP - 1) are annotated in dbSNP database." > ${out}/all_variants_common_summary.txt

echo "===Step 4: Annotating based on ClinVar Database==="
#Annotate using SNPSift with reference to CLINVAR Database
java -Xmx20g -jar ~/Softwares/snpEff/SnpSift.jar annotate ~/Softwares/snpEff/data/clinvar/clinvar_20250321.vcf.gz ${out}all_variants.common.vcf > ${out}all_variants.common.clinvar.vcf
bgzip -c ${out}all_variants.common.clinvar.vcf > ${out}all_variants.common.clinvar.vcf.gz
rm ${out}all_variants.common.clinvar.vcf

bcftools view -i 'INFO/CLNSIG="Pathogenic"' ${out}all_variants.common.clinvar.vcf.gz -Oz | bcftools query -f '[%CHROM %POS %TYPE %ID %GT %GQ %INFO/CLNDN %INFO/CLNSIG %INFO/CLNSIGCONF %INFO/CLNSIGINCL %INFO/CLNVC %INFO/GENEINFO\n]' > ${out}all_variants.common.clinvar_patho.txt

bcftools view -i 'INFO/CLNSIG="Likely_pathogenic"' ${out}all_variants.common.clinvar.vcf.gz -Oz | bcftools query -f '[%CHROM %POS %TYPE %ID %GT %GQ %INFO/CLNDN %INFO/CLNSIG %INFO/CLNSIGCONF %INFO/CLNSIGINCL %INFO/CLNVC %INFO/GENEINFO\n]' > ${out}all_variants.common.clinvar_likelypatho.txt

bcftools view -i 'INFO/ONC="Oncogenic"' ${out}all_variants.common.clinvar.vcf.gz -Oz | bcftools query -f '[%CHROM %POS %TYPE %ID %GT %GQ %INFO/CLNDN %INFO/CLNSIG %INFO/CLNSIGCONF %INFO/CLNSIGINCL %INFO/CLNVC %INFO/GENEINFO\n]' > ${out}all_variants.common.clinvar_oncogenic.txt

count_patho=$(awk -F ' ' '{print $7}' ${out}all_variants.common.clinvar_patho.txt | sort -u | wc -l)
count_likelypatho=$(awk -F ' ' '{print $7}' ${out}all_variants.common.clinvar_likelypatho.txt | sort -u | wc -l)
count_oncogenic=$(awk -F ' ' '{print $7}' ${out}all_variants.common.clinvar_oncogenic.txt | sort -u | wc -l)

echo "The input callset has $count_patho pathogenic variants, $count_likelypatho likely pathogenic variants, and $oncogenic oncogenic variants based on ClinVar database." >> ${out}all_variants_common_summary.txt

echo "Completed"
