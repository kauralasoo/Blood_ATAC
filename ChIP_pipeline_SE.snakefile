import uuid
import os

#Align reads to the reference genome using BWA-ALN (<70bp SE reads)
rule align_with_bwa_aln:
	input:
		"processed/{dataset}/fastq/{sample}.fastq.gz"
	output:
		"processed/{dataset}/aligned/{sample}.bam"
	params:
		genome = config["bwa_index"],
		rg="@RG\tID:{sample}\tSM:{sample}",
		tmp_fq = "/tmp/" + uuid.uuid4().hex + ".fastq.gz",
		tmp_sai = "/tmp/" + uuid.uuid4().hex + ".sai",
		tmp_bam = "/tmp/" + uuid.uuid4().hex + ".bam"
	resources:
		mem = 12000
	threads: 4
	shell:
		"""
		module load samtools-1.6
		module load bwa-0.7.12
		cp {input} {params.tmp_fq1}
		bwa aln -t {threads} {params.genome} {params.tmp_fq} > {params.tmp_sai}
		bwa samse -r '{params.rg}' {params.genome} {params.tmp_sai} {params.tmp_fq} | samtools view -b - > {params.tmp_bam}
		cp {params.tmp_bam} {output}
		rm {params.tmp_bam}
		rm {params.tmp_sai}
		rm {params.tmp_fq}
		"""

#Make sure that all final output files get created
rule make_all:
	input:
		expand("processed/{{dataset}}/aligned/{sample}.bam", sample=config["samples"])
	output:
		"processed/{dataset}/out.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done' > {output}"

