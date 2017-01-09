configfile: "config.yaml"

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
		"python scripts/extractNexteraBarcode.py --fastq {input.f1} --type read1 > {output.f1} && "
		"python scripts/extractNexteraBarcode.py --fastq {input.f2} --type read2 > {output.f2}"

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
		genome = config["bwa_index"],
		rg="@RG\tID:{sample}\tSM:{sample}"
	resources:
		mem = 12000
	threads: 4
	shell:
		"bwa mem -M -R '{params.rg}' -t {threads} {params.genome} {input} | samtools view -b - > {output}"

#Sort BAM files by coordinates
rule sort_bams_by_position:
	input:
		"processed/aligned/{sample}.bam"
	output:
		"processed/aligned/{sample}.sorted.bam"
	resources:
		mem = 4000
	threads: 4
	shell:
		"samtools sort -T processed/aligned/{wildcards.sample} -O bam -@ {threads} {input} > {output}"

#Index sorted bams
rule index_bams:
	input:
		"processed/aligned/{sample}.sorted.bam"
	output:
		"processed/aligned/{sample}.sorted.bam.bai"
	resources:
		mem = 50
	threads: 1
	shell:
		"samtools index {input}"

#Keep only properly paired reads from nuclear chromosomes (remove MT)
rule filter_properly_paired:
	input:
		bam = "processed/aligned/{sample}.sorted.bam",
		index = "processed/aligned/{sample}.sorted.bam.bai"
	output:
		"processed/filtered/{sample}.filtered.bam"
	resources:
		mem = 100
	threads: 1
	params:
		chr_list = "1 10 11 12 13 14 15 16 17 18 19 2 20 21 22 3 4 5 6 7 8 9 X Y"
	shell:
		"samtools view -h -b -f 2 {input} {params.chr_list} > {output}"

#Remove BWA entry from the BAM file header (conflicts with MarkDuplicates)
rule remove_bwa_header:
	input:
		"processed/filtered/{sample}.filtered.bam"
	output:
		new_header = "processed/filtered/{sample}.new_header.txt",
		bam = "processed/filtered/{sample}.reheadered.bam"
	resources:
		mem = 100
	threads: 1
	shell:
		"samtools view -H {input} | grep -v 'ID:bwa' > {output.new_header} &&"
		"samtools reheader {output.new_header} {input} > {output.bam}"

#Remove duplicates using Picard MarkDuplicates
rule remove_duplicates:
	input:
		"processed/filtered/{sample}.reheadered.bam"
	output:
		bam = "processed/filtered/{sample}.no_duplicates.bam",
		metrics = "processed/metrics/{sample}.MarkDuplicates.txt"
	resources:
		mem = 2200
	threads: 4
	shell:
		"{config[picard_path]}r MarkDuplicates I={input} O={output.bam} REMOVE_DUPLICATES=true METRICS_FILE= {output.metrics}"

#Count the number of reads per chromosome (QC metric)
rule reads_per_chromosome:
	input:
		"processed/aligned/{sample}.sorted.bam"
	output:
		"processed/metrics/{sample}.chr_counts.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"samtools view {input} | cut -f3 | sort | uniq -c > {output}"


#Sort BAM files by coordinates
rule sort_bams_by_name:
	input:
		"processed/filtered/{sample}.no_duplicates.bam"
	output:
		"processed/filtered/{sample}.sortedByName.bam"
	resources:
		mem = 4000
	threads: 4
	shell:
		"samtools sort -T processed/filtered/{wildcards.sample} -O bam -@ {threads} -n {input} > {output}"

#Convert the bam file to a bed file of fragments
#First to bedpe and then to bed
rule bam_to_fragment_bed:
	input:
		"processed/filtered/{sample}.sortedByName.bam"
	output:
		"processed/filtered/{sample}.fragments.bed.gz"
	resources:
		mem = 1000
	threads: 2
	shell:
		"bedtools bamtobed -bedpe -i {input} | python scripts/bedpe2bed.py --maxFragmentLength 1000 | sort -k 1,1 | gzip > {output}"

#Count the number of occurences of each fragment length in the fragments bed file
rule count_fragment_lengths:
	input:
		"processed/filtered/{sample}.fragments.bed.gz"
	output:
		"processed/metrics/{sample}.fragment_lengths.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"zcat {input} | cut -f5 | sort -n | uniq -c > {output}"

#Convert fragment bed file into bigwig
rule convert_bed_to_bigwig:
	input:
		"processed/filtered/{sample}.fragments.bed.gz"
	output:
		bedgraph = "processed/bigwig/{sample}.bg.gz",
		bigwig = "processed/bigwig/{sample}.bw"
	params:
		bedgraph = "processed/bigwig/{sample}.bg"
	resources:
		mem = 3000
	threads: 1
	shell:
		"bedtools genomecov -bga -i {input} -g {config[chromosome_lengths]} > {params.bedgraph} && "
		"bedGraphToBigWig {params.bedgraph} {config[chromosome_lengths]} {output.bigwig} && "
		"gzip {params.bedgraph}"


rule make_all_bigwigs:
	input:
		expand("processed/bigwig/{sample}.bw", sample=config["samples"])
	output:
		"processed/out.txt"
	shell:
		"echo 'Done' > {output}"

