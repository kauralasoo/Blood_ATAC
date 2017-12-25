import uuid
import os

#Align reads to the reference genome using BWA-ALN (<70bp SE reads)
rule align_with_bwa_aln:
	input:
		"processed/{dataset}/fastq/{sample}.fastq.gz"
	output:
		temp("processed/{dataset}/aligned/{sample}.bam")
	params:
		index_dir = config["bwa_index_dir"],
		index_name = config["bwa_index_name"],
		rg="@RG\tID:{sample}\tSM:{sample}",
		tmp_fq = "/tmp/" + uuid.uuid4().hex + ".fastq.gz",
		tmp_sai = "/tmp/" + uuid.uuid4().hex + ".sai",
		tmp_bam = "/tmp/" + uuid.uuid4().hex + ".bam",
		tmp_index_dir = "/tmp/" + uuid.uuid4().hex + "/"
	resources:
		mem = 12000
	threads: 8
	shell:
		"""
		module load samtools-1.6
		module load bwa-0.7.12
		cp {input} {params.tmp_fq}
		mkdir {params.tmp_index_dir}
		cp {params.index_dir}/* {params.tmp_index_dir}
		bwa aln -t {threads} {params.tmp_index_dir}{params.index_name} {params.tmp_fq} > {params.tmp_sai}
		bwa samse -r '{params.rg}' {params.tmp_index_dir}{params.index_name} {params.tmp_sai} {params.tmp_fq} | samtools view -b - > {params.tmp_bam}
		cp {params.tmp_bam} {output}
		rm {params.tmp_bam}
		rm {params.tmp_sai}
		rm {params.tmp_fq}
		rm -r {params.tmp_index_dir}
		"""

#Sort BAM files by coordinates
rule sort_bams_by_position:
	input:
		"processed/{dataset}/aligned/{sample}.bam"
	output:
		protected("processed/{dataset}/sorted_bam/{sample}.sortedByCoords.bam")
	params:
		local_tmp = "/tmp/" + uuid.uuid4().hex + "/"
	resources:
		mem = 8000
	threads: 4
	shell:
		"""
		module load samtools-1.6
		mkdir {params.local_tmp}
		cp {input} {params.local_tmp}/{wildcards.sample}.bam
		samtools sort -T {params.local_tmp}/{wildcards.sample} -O bam -@ {threads} -o {params.local_tmp}/{wildcards.sample}.sortedByCoords.bam {params.local_tmp}/{wildcards.sample}.bam
		cp {params.local_tmp}/{wildcards.sample}.sortedByCoords.bam {output}
		rm -r {params.local_tmp}
		"""

#Index sorted bams
rule index_bams:
	input:
		"processed/{dataset}/sorted_bam/{sample}.sortedByCoords.bam"
	output:
		"processed/{dataset}/sorted_bam/{sample}.sortedByCoords.bam.bai"
	resources:
		mem = 50
	threads: 1
	shell:
		"""
		module load samtools-1.6
		samtools index {input}
		"""

#Keep only properly paired reads from nuclear chromosomes (remove MT)
rule filter_properly_paired:
	input:
		bam = "processed/{dataset}/sorted_bam/{sample}.sortedByCoords.bam",
		index = "processed/{dataset}/sorted_bam/{sample}.sortedByCoords.bam.bai"
	output:
		temp("processed/{dataset}/filtered/{sample}.filtered.bam")
	resources:
		mem = 100
	threads: 1
	params:
		chr_list = "1 10 11 12 13 14 15 16 17 18 19 2 20 21 22 3 4 5 6 7 8 9 X Y"
	shell:
		"""
		module load samtools-1.6
		samtools view -h -b {input.bam} {params.chr_list} > {output}
		"""

#Remove BWA entry from the BAM file header (conflicts with MarkDuplicates)
rule remove_bwa_header:
	input:
		"processed/{dataset}/filtered/{sample}.filtered.bam"
	output:
		new_header = temp("processed/{dataset}/filtered/{sample}.new_header.txt"),
		bam = temp("processed/{dataset}/filtered/{sample}.reheadered.bam")
	resources:
		mem = 100
	threads: 1
	shell:
		"""
		module load samtools-1.6
		samtools view -H {input} | grep -v 'ID:bwa' > {output.new_header}
		samtools reheader {output.new_header} {input} > {output.bam}
		"""

#Remove duplicates using Picard MarkDuplicates
rule remove_duplicates:
	input:
		"processed/{dataset}/filtered/{sample}.reheadered.bam"
	output:
		bam = protected("processed/{dataset}/filtered/{sample}.no_duplicates.bam"),
		metrics = protected("processed/{dataset}/metrics/{sample}.MarkDuplicates.txt")
	resources:
		mem = 2200
	threads: 4
	shell:
		"""
		module load jdk-1.8.0_25
		java -jar {config[picard_path]} MarkDuplicates I={input} O={output.bam} REMOVE_DUPLICATES=true METRICS_FILE={output.metrics}
		"""


#Make sure that all final output files get created
rule make_all:
	input:
		expand("processed/{{dataset}}/aligned/{sample}.bam", sample=config["samples"]),
		expand("processed/{{dataset}}/sorted_bam/{sample}.sortedByCoords.bam.bai", sample=config["samples"]),
		expand("processed/{{dataset}}/filtered/{sample}.no_duplicates.bam", sample=config["samples"])
	output:
		"processed/{dataset}/out.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done' > {output}"

