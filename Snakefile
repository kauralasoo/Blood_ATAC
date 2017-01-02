SAMPLES = ["HSC_4983_1"]

rule SRA_to_fastq:
	input:
		"processed/SRR/{sample}.sra"
	output:
		"processed/fastq/{sample}_1.fastq.gz",
		"processed/fastq/{sample}_2.fastq.gz"
	shell:
		"fastq-dump --split-files --gzip --skip-technical --readids --dumpbase --clip --outdir processed/fastq/ {input}"

rule rename_fastq:
	input:
		f1 = "processed/fastq/{sample}_1.fastq.gz",
		f2 = "processed/fastq/{sample}_2.fastq.gz"
	output:
		f1 = "processed/fastq/{sample}.1.fastq.gz",
		f2 = "processed/fastq/{sample}.2.fastq.gz"
	shell:
		"mv {input.f1} {output.f1} && mv {input.f2} {output.f2}" 
