
#Concat vcf files
zcat full_header51.txt.gz chr1.gz chr11.gz chr12.gz chr13.gz chr14.gz chr15.gz chr16.gz chr17.gz chr18.gz chr19.gz chr2.gz chr20.gz chr21.gz chr22.gz | bgzip > CTCF_51_samples.GRCh37.vcf.gz

#Update GT and DS fields in the header
reheader -h new_header.txt CTCF_51_samples.GRCh37.vcf.gz > CTCF_51_samples.GRCh37.reheadered.vcf.gz

#Run Snakemake
snakemake --cluster ../../scripts/snakemake_submit_UT.py -np --snakefile process_genotypes.snakefile


snakemake --cluster ../../scripts/snakemake_submit_UT.py -np --snakefile process_genotypes.snakefile CTCF/genotypes/vcf/GRCh38/CTCF_51_samples.GRCh38.common.sorted.vcf.gz --jobs 1


#Keep only common SNPs
bcftools filter -i 'MAF[0] >= 0.05' -O z CTCF_51_samples.GRCh38.vcf.gz > CTCF_51_samples.GRCh38.common.vcf.gz