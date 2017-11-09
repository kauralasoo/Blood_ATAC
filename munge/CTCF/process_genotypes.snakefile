rule rename_chromosomes:
    input:
        vcf = "CTCF/genotypes/vcf/GRCh37/CTCF_51_samples.GRCh37.vcf.gz",
        chromosome_map = "../../data/liftOver/GRCh38ToHg38_chromosome_map.txt",
    output:
        vcf = "CTCF/genotypes/vcf/GRCh38/CTCF_51_samples.hg19.vcf.gz",
    threads: 1
    resources:
        mem = 1000
    shell:
        """
        module load bcftools-1.6
        bcftools annotate -O z --rename-chrs {input.chromosome_map} {input.vcf} > {output.vcf}
        """
