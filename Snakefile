SAMPLES = ["HSC_4983_1"]

#Convert SRA files to fastq
rule SRA_to_fastq:
	input:
		"processed/SRR/{sample}.sra"
	output:
		"processed/fastq/{sample}_1.fastq.gz",
		"processed/fastq/{sample}_2.fastq.gz"
	resources:
		mem = 1000
	threads: 1
	shell:
		"fastq-dump --split-files --gzip --skip-technical --dumpbase --clip --outdir processed/fastq/ {input}"

#Rename the fastq files to specify pairing correctly
rule rename_fastq:
	input:
		f1 = "processed/fastq/{sample}_1.fastq.gz",
		f2 = "processed/fastq/{sample}_2.fastq.gz"
	output:
		f1 = "processed/fastq/{sample}.1.fastq.gz",
		f2 = "processed/fastq/{sample}.2.fastq.gz"
	resources:
		mem = 50
	threads: 1
	shell:
		"mv {input.f1} {output.f1} && mv {input.f2} {output.f2}" 

#Extract Nextera barcodes from the fastq files
rule extract_nextera_barcodes:
	input:
		f1 = "processed/fastq/{sample}.1.fastq.gz",
		f2 = "processed/fastq/{sample}.2.fastq.gz"
	output:
		f1 = "processed/trimmed/{sample}.1.barcode.txt",
		f2 = "processed/trimmed/{sample}.2.barcode.txt"
	resources:
		mem = 50
	threads: 1
	shell:
		"python ~/software/utils/fastq/extractNexteraBarcode.py --fastq {input.f1} --type read1 > {output.f1} && "
		"python ~/software/utils/fastq/extractNexteraBarcode.py --fastq {input.f2} --type read2 > {output.f2}"

#Trim adapters from fastq files
rule trim_adapters:
	input:
		fastq1 = "processed/fastq/{sample}.1.fastq.gz",
		fastq2 = "processed/fastq/{sample}.2.fastq.gz",
		barcode1 = "processed/trimmed/{sample}.1.barcode.txt",
		barcode2 = "processed/trimmed/{sample}.2.barcode.txt"
	params:
		prefix="processed/trimmed/{sample}"
	resources:
		mem = 1000
	threads: 1
	output:
		"processed/trimmed/{sample}-trimmed-pair1.fastq.gz",
		"processed/trimmed/{sample}-trimmed-pair2.fastq.gz"
	script:
		"scripts/trim_adapters.py"

#Rename trimmed fastq files from skewer
rule rename_trimmed_fastq:
	input:
		f1 = "processed/trimmed/{sample}-trimmed-pair1.fastq.gz",
		f2 = "processed/trimmed/{sample}-trimmed-pair2.fastq.gz"
	output:
		f1 = "processed/trimmed/{sample}.1.trimmed.fastq.gz",
		f2 = "processed/trimmed/{sample}.2.trimmed.fastq.gz"
	resources:
		mem = 50
	threads: 1
	shell:
		"mv {input.f1} {output.f1} && mv {input.f2} {output.f2}" 

#Align reads to the reference genome using BWA
rule align_reads:
	input:
		"processed/trimmed/{sample}.1.trimmed.fastq.gz",
		"processed/trimmed/{sample}.2.trimmed.fastq.gz"
	output:
		"processed/aligned/{sample}.bam"
	params:
		genome = "../../../annotations/GRCh38/bwa_index/GRCh38",
		rg="@RG\tID:{sample}\tSM:{sample}"
	resources:
		mem = 12000
	threads: 4
	shell:
		"bwa mem -M -R '{params.rg}' -t {threads} {params.genome} {input} | samtools view -b - > {output}"



