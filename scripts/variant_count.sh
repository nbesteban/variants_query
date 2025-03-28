#!/bin/sh

#$1 and $2 should be passed when calling the script. $1 is the vcf file while $2 is the output file
chromlist=$(cat data/chromosomes.txt)

for chrom in $chromlist[@];
do
	echo "Counting all variants of chromosome... $chrom"
	count_all=$(bcftools view -r $chrom $1 | grep -v -c '^#')
	count_novel=$(bcftools view -r $chrom -n $1 | grep -v -c '^#')
	echo "Counting bi-allelic SNP variants of chromosome... $chrom"
	count_bisnps=$(bcftools view -r $chrom -m2 -M2 -v snps $1 | grep -v -c '^#')
	count_novel_bisnps=$(bcftools view -r $chrom -m2 -M2 -v snps -n $1 | grep -v -c '^#')
	echo "Counting tri-allelic SNP variants of chromosome... $chrom"
	count_trisnps=$(bcftools view -r $chrom -m3 -M3 -v snps $1 | grep -v -c '^#')
	count_novel_trisnps=$(bcftools view -r $chrom -m3 -M3 -v snps -n $1 | grep -v -c '^#')
	echo "Counting Quad-allelic SNP variants of chromosome $chrom..."
	count_quadsnps=$(bcftools view -r $chrom -m4 -M4 -v snps $1 | grep -v -c '^#')
        count_novel_quadsnps=$(bcftools view -r $chrom -m4 -M4 -v snps -n $1 | grep -v -c '^#')
	echo "Counting INDEL variants of chromosome $chrom..."
	count_in=$(bcftools view -r $chrom --types indels --include 'ILEN>0' $1 | grep -v -c '^#')
	count_novel_in=$(bcftools view -r $chrom --types indels -n --include 'ILEN>0' $1 | grep -v -c '^#')
	count_del=$(bcftools view -r $chrom --types indels --include 'ILEN<0' $1 | grep -v -c '^#')
	count_novel_del=$(bcftools view -r $chrom --types indels -n --include 'ILEN<0' $1 | grep -v -c '^#')
	echo "Printing summary of the variant counts for chromosome $chrom."
	echo -e "${chrom}\t${count_all}\t${count_novel}\t${count_bisnps}\t${count_novel_bisnps}\t${count_trisnps}\t${count_novel_trisnps}\t${count_quadsnps}\t${count_novel_quadsnps}\t${count_in}\t${count_novel_in}\t${count_del}\t${count_novel_del}" >> $2count_perchrom.txt
done

sed -i $'1iChr\tTotal_Variants\tNovelSNPs\tBiAlllelic_SNPs\tNovel_BiAllelic_SNPs\tTriAlllelic_SNPs\tNovel_TriAllelic_SNPs\tQuadAlllelic_SNPs\tNovel_QuadAllelic_SNPs\tInsertion\tNovel_Insertion\tDeletion\tNover_Deletion' $2count_perchrom.txt
SUM=$(awk '{for (i=1;i<=NF;i++) sum[i]+=$i;}; END{for (i in sum) print sum[i]}' ORS="\t" ${2}count_perchrom.txt)
