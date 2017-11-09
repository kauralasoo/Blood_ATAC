rule rename_chromosomes:
    input:
        vcf = "CTCF/genotypes/vcf/GRCh37/CTCF_51_samples.GRCh37.reheadered.vcf.gz",
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

rule run_CrossMap:
    input:
        vcf = "CTCF/genotypes/vcf/GRCh38/CTCF_51_samples.hg19.vcf.gz",
        chain = "../../data/liftOver/hg19ToHg38.over.chain",
        ref_genome = "/gpfs/rocket/home/a72094/annotations/hg38/hg38.fa",
    output:
        vcf = "CTCF/genotypes/vcf/GRCh38/CTCF_51_samples.hg38.vcf.gz",
    params:
        temp_vcf = "CTCF/genotypes/vcf/GRCh38/CTCF_51_samples.hg38.vcf"
    threads: 1
    resources:
        mem = 1000
    shell:
        """
        source activate py2.7
        CrossMap.py vcf {input.chain} {input.vcf} {input.ref_genome} {params.temp_vcf}
        bgzip {params.temp_vcf}
        """

rule postprocess_CrossMap:
    input:
        vcf = "CTCF/genotypes/vcf/GRCh38/CTCF_51_samples.hg38.vcf.gz"
    output:
        vcf = "CTCF/genotypes/vcf/GRCh38/CTCF_51_samples.hg38.post.vcf.gz"
    threads: 1
    resources:
        mem = 1000
    shell:
        """
        source activate py2.7
        module load bcftools-1.6
        python ../../scripts/postprocessCrossmap.py --vcf {input.vcf} | bgzip > {output.vcf}
        """