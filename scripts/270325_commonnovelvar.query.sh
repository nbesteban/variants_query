#!/bin/sh

# Ensure correct argument usage
if [ "$1" != "--raw" ] || [ "$3" != "--out" ]; then
	    echo "Usage: $0 --raw <input.vcf.gz> --out <OutputDirectory>"
	        exit 1
fi

# Command-line input arguments
raw="$2"	#Path to the raw vcf (compressed) dataset
out="$4"	#Path to the directory for the output file

echo "Investigating Novel Variants"
echo "====Step 1: Extracting novel variants from the data==="
#Extract novel variants from the data, those that have no rsid in the dbSNPs (annotated by '.' in chr position).
#bcftools view -n $raw -Oz -o ${out}novel_variants.vcf.gz
#tabix -p vcf ${out}novel_variants.vcf.gz

echo "===Step 2: Removing related individuals, low quality variants, and samples with high-missing data==="
bcftools view -S ^samplestoremove.list -i 'MAF>0.01 && QUAL>30 && FMT/GQ>10 && FILTER="PASS"' ${out}novel_variants.vcf.gz -Oz -o ${out}novel_variants.common.vcf.gz
tabix -p vcf ${out}novel_variants.common.vcf.gz
bgzip -d ${out}novel_variants.common.vcf.gz > ${out}novel_variants.common.vcf

echo "===Step 3: Annotating using SnpEff==="
#Annotate the common novel variants (MAF>0.01) using SnpEff
java -Xmx20g -jar ~/Softwares/snpEff/snpEff.jar -v hg38 ${out}novel_variants.common.vcf.gz -stats ${out}novel_variants.common.snpEff > ${out}novel_variants.common.snpEff.vcf
bgzip -c ${out}novel_variants.common.snpEff.vcf > ${out}novel_variants.common.snpEff.vcf.gz
rm ${out}novel_variants.common.snpEff.vcf

echo "===Step 4: Annotating using SIFT==="
#Annotate using SIFT, single transcipt only
mkdir ${out}novel_variants.common.sift
java -jar ~/Softwares/SIFT4G_Annotator.jar -c -i ${out}novel_variants.common.vcf -d ../resources/GRCh38.83.chr/ -r ${out}novel_variants.common.sift
bgzip -c ${out}novel_variants.common.sift/novel_variants.common_SIFTpredictions.vcf > ${out}novel_variants.common.sift/novel_variants.common_SIFTpredictions.vcf.gz
rm ${out}novel_variants.common.sift/novel_variants.common_SIFTpredictions.vcf

sed -n '1p;/DELETERIOUS/p' ${out}novel_variants.common.sift/novel_variants.common_SIFTannotations.xls | grep -v "*WARNING!" > ${out}novel_variants.common.sift_deleterious.txt
sed -n '1p;/novel/p' ${out}novel_variants.common.sift_deleterious.txt > ${out}novel_variants.common.sift_deleterious_novel.txt
sed -n '1p;/novel/!p' ${out}novel_variants.common.sift_deleterious.txt > ${out}novel_variants.common.sift_deleterious_dbSNP.txt

var_count_del=$(cat ${out}novel_variants.common.sift_deleterious.txt | wc -l)
var_count_del_novel=$(cat ${out}novel_variants.common.sift_deleterious_novel.txt | wc -l)
var_count_del_dbSNP=$(cat ${out}novel_variants.common.sift_deleterious_dbSNP.txt | wc -l)
count_del=$(awk -F '\t' '{print $6}' ${out}novel_variants.common.sift_deleterious.txt | sort -u | wc -l)
count_del_novel=$(awk -F '\t' '{print $6}' ${out}novel_variants.common.sift_deleterious_novel.txt | sort -u | wc -l)
count_del_dbSNP=$(awk -F '\t' '{print $6}' ${out}novel_variants.common.sift_deleterious_dbSNP.txt | sort -u | wc -l)

#need to subtract 1 to remove the header line count since we counted it including the header line
echo "The input callset has $(expr $var_count_del - 1) deleterious variants that transcribed to $(expr $count_del - 1) genes. These deleterious variants consists of $(expr $var_count_del_novel - 1) novel variants to the dataset and $(expr $var_count_del_dbSNP - 1) already annotated in dbSNP database, which transcribed to $(expr $count_del_novel - 1) and $(expr $count_del_dbSNP - 1) genes, respectively." > ${out}/novel_variants_summary.txt

echo "===Step 5: Annotating based on ClinVar Database==="
#Annotate using SNPSift with reference to CLINVAR Database
java -Xmx20g -jar ~/Softwares/snpEff/SnpSift.jar annotate ~/Softwares/snpEff/data/clinvar/clinvar_20250321.vcf.gz ${out}novel_variants.common.vcf > ${out}novel_variants.common.clinvar.vcf
bgzip -c ${out}novel_variants.common.clinvar.vcf > ${out}novel_variants.common.clinvar.vcf.gz
rm ${out}novel_variants.common.clinvar.vcf

bcftools view -i 'INFO/CLNSIG="Pathogenic"' ${out}novel_variants.common.clinvar.vcf.gz -Oz | bcftools query -f '[%CHROM %POS %TYPE %ID %GT %GQ %INFO/CLNDN %INFO/CLNSIG %INFO/CLNSIGCONF %INFO/CLNSIGINCL %INFO/CLNVC %INFO/GENEINFO %AF\n]' > ${out}novel_variants.common.clinvar_patho.txt

bcftools view -i 'INFO/CLNSIG="Likely_pathogenic"' ${out}novel_variants.common.clinvar.vcf.gz -Oz | bcftools query -f '[%CHROM %POS %TYPE %ID %GT %GQ %INFO/CLNDN %INFO/CLNSIG %INFO/CLNSIGCONF %INFO/CLNSIGINCL %INFO/CLNVC %INFO/GENEINFO %AF\n]' > ${out}novel_variants.common.clinvar_likelypatho.txt

bcftools view -i 'INFO/CLNSIG="Benign"' ${out}novel_variants.common.clinvar.vcf.gz -Oz | bcftools query -f '[%CHROM %POS %TYPE %ID %GT %GQ %INFO/CLNDN %INFO/CLNSIG %INFO/CLNSIGCONF %INFO/CLNSIGINCL %INFO/CLNVC %INFO/GENEINFO %AF\n]' > ${out}novel_variants.common.clinvar_benign.txt

bcftools view -i 'INFO/CLNSIG="Likely_benign"' ${out}novel_variants.common.clinvar.vcf.gz -Oz | bcftools query -f '[%CHROM %POS %TYPE %ID %GT %GQ %INFO/CLNDN %INFO/CLNSIG %INFO/CLNSIGCONF %INFO/CLNSIGINCL %INFO/CLNVC %INFO/GENEINFO %AF\n]' > ${out}novel_variants.common.clinvar_likelybenign.txt

bcftools view -i 'INFO/ONC="Oncogenic"' ${out}novel_variants.common.clinvar.vcf.gz -Oz | bcftools query -f '[%CHROM %POS %TYPE %ID %GT %GQ %INFO/CLNDN %INFO/CLNSIG %INFO/CLNSIGCONF %INFO/CLNSIGINCL %INFO/CLNVC %INFO/GENEINFO %AF\n]' > ${out}novel_variants.common.clinvar_oncogenic.txt



count_patho=$(awk -F ' ' '{print $1}' ${out}novel_variants.common.clinvar_patho.txt | sort -u | wc -l)
count_likelypatho=$(awk -F ' ' '{print $1}' ${out}novel_variants.common.clinvar_likelypatho.txt | sort -u | wc -l)
count_benign=$(awk -F ' ' '{print $1}' ${out}novel_variants.common.clinvar_benign.txt | sort -u | wc -l)
count_likelybenign=$(awk -F ' ' '{print $1}' ${out}novel_variants.common.clinvar_likelybenign.txt | sort -u | wc -l)
count_oncogenic=$(awk -F ' ' '{print $1}' ${out}novel_variants.common.clinvar_oncogenic.txt | sort -u | wc -l)


echo "The input callset has $count_patho pathogenic variants, $count_likelypatho likely pathogenic variants, $count_benign benign variants, $count_likelybenign likely benign variants, and $oncogenic oncogenic variants based on ClinVar database." >> ${out}novel_variants_summary.txt

echo "Completed"
