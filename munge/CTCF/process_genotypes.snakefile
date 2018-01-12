rule rename_chromosomes:
    input:
        vcf = "genotypes/vcf/GRCh37/{study}.GRCh37.vcf.gz",
        chromosome_map = "data/liftOver/GRCh38ToHg38_chromosome_map.txt",
    output:
        vcf = "genotypes/vcf/GRCh38/{study}.hg19.vcf.gz",
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
        vcf = "genotypes/vcf/GRCh38/{study}.hg19.vcf.gz",
        chain = "data/liftOver/hg19ToHg38.over.chain",
        ref_genome = "/gpfs/rocket/home/a72094/annotations/hg38/hg38.fa",
    output:
        vcf = "genotypes/vcf/GRCh38/{study}.hg38.vcf.gz",
    params:
        temp_vcf = "genotypes/vcf/GRCh38/{study}.hg38.vcf"
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
        vcf = "genotypes/vcf/GRCh38/{study}.hg38.vcf.gz"
    output:
        vcf = "genotypes/vcf/GRCh38/{study}.hg38.post.vcf.gz"
    threads: 1
    resources:
        mem = 1000
    shell:
        """
        source activate py2.7
        module load bcftools-1.6
        python scripts/postprocessCrossmap.py --vcf {input.vcf} | bgzip > {output.vcf}
        """

rule rename_chromosomes_back:
    input:
        vcf = "genotypes/vcf/GRCh38/{study}.hg38.post.vcf.gz",
        chromosome_map = "data/liftOver/Hg38ToGRCh38_chromosome_map.txt",
    output:
        vcf = "genotypes/vcf/GRCh38/{study}.GRCh38.vcf.gz"
    threads: 1
    resources:
        mem = 1000
    shell:
        """
        module load bcftools-1.6
        bcftools annotate -O z --rename-chrs {input.chromosome_map} {input.vcf} > {output.vcf}
        """

rule sort_vcf:
    input:
        vcf = "genotypes/vcf/GRCh38/{study}.GRCh38.vcf.gz"
    output:
        vcf = "genotypes/vcf/GRCh38/{study}.GRCh38.sorted.vcf.gz"
    threads: 1
    resources:
        mem = 35000
    shell:
        """
        module load bcftools-1.6
        bcftools sort -m 20000M -o {output.vcf} -O z {input.vcf}
        """

rule filter_ref_allele:
    input:
        vcf = "genotypes/vcf/GRCh38/{study}.GRCh38.sorted.vcf.gz",
        fasta = "/gpfs/rocket/home/a72094/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    output:
        vcf = "genotypes/vcf/GRCh38/{study}.GRCh38.sorted.ref.vcf.gz"
    threads: 1
    resources:
        mem = 3000
    shell:
        """
        module load bcftools-1.6
        bcftools norm -c x -O z -f {input.fasta} {input.vcf} > {output.vcf}
        """

rule remove_multialleic:
    input:
        vcf = "genotypes/vcf/GRCh38/{study}.GRCh38.sorted.ref.vcf.gz"
    output:
        vcf = "genotypes/vcf/GRCh38/{study}.GRCh38.sorted.ref.filtered.vcf.gz"
    threads: 1
    resources:
        mem = 2000
    shell:
        """
        module load bcftools-1.6
        bcftools norm -m+any {input.vcf} | bcftools view -m2 -M2 - | bcftools annotate --set-id +'%CHROM\_%POS' | bcftools norm -d both -O z > {output.vcf}
        """

rule keep_common:
    input:
        "genotypes/vcf/GRCh38/{study}.GRCh38.sorted.ref.filtered.vcf.gz"
    output:
        "genotypes/vcf/GRCh38/{study}.GRCh38.final.vcf.gz"
    threads: 1
    resources:
        mem = 3000
    shell:
        """
        module load bcftools-1.6
        bcftools filter -i 'MAF[0] >= 0.05' -O z {input} > {output}
        """
