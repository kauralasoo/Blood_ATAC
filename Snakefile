configfile: "config.yaml"

#Convert SRA files to fastq
rule SRA_to_fastq:
	input:
		"{dataset}/SRR/{sample}.sra"
	output:
		"{dataset}/fastq/{sample}_1.fastq.gz",
		"{dataset}/fastq/{sample}_2.fastq.gz"
	resources:
		mem = 1000
	threads: 1
	shell:
		"fastq-dump --split-files --gzip --skip-technical --dumpbase --clip --outdir {wildcards.dataset}/fastq/ {input}"

#Rename the fastq files to specify pairing correctly
rule rename_fastq:
	input:
		f1 = "{dataset}/fastq/{sample}_1.fastq.gz",
		f2 = "{dataset}/fastq/{sample}_2.fastq.gz"
	output:
		f1 = "{dataset}/fastq/{sample}.1.fastq.gz",
		f2 = "{dataset}/fastq/{sample}.2.fastq.gz"
	resources:
		mem = 50
	threads: 1
	shell:
		"mv {input.f1} {output.f1} && mv {input.f2} {output.f2}" 

#Extract Nextera barcodes from the fastq files
rule extract_nextera_barcodes:
	input:
		f1 = "{dataset}/fastq/{sample}.1.fastq.gz",
		f2 = "{dataset}/fastq/{sample}.2.fastq.gz"
	output:
		f1 = "{dataset}/trimmed/{sample}.1.barcode.txt",
		f2 = "{dataset}/trimmed/{sample}.2.barcode.txt"
	resources:
		mem = 50
	threads: 1
	shell:
		"python scripts/extractNexteraBarcode.py --fastq {input.f1} --type read1 > {output.f1} && "
		"python scripts/extractNexteraBarcode.py --fastq {input.f2} --type read2 > {output.f2}"

#Trim adapters from fastq files
rule trim_adapters:
	input:
		fastq1 = "{dataset}/fastq/{sample}.1.fastq.gz",
		fastq2 = "{dataset}/fastq/{sample}.2.fastq.gz",
		barcode1 = "{dataset}/trimmed/{sample}.1.barcode.txt",
		barcode2 = "{dataset}/trimmed/{sample}.2.barcode.txt"
	params:
		prefix="{dataset}/trimmed/{sample}"
	resources:
		mem = 1000
	threads: 3
	output:
		"{dataset}/trimmed/{sample}-trimmed-pair1.fastq.gz",
		"{dataset}/trimmed/{sample}-trimmed-pair2.fastq.gz"
	script:
		"scripts/trim_adapters.py"

#Rename trimmed fastq files from skewer
rule rename_trimmed_fastq:
	input:
		f1 = "{dataset}/trimmed/{sample}-trimmed-pair1.fastq.gz",
		f2 = "{dataset}/trimmed/{sample}-trimmed-pair2.fastq.gz"
	output:
		f1 = "{dataset}/trimmed/{sample}.1.trimmed.fastq.gz",
		f2 = "{dataset}/trimmed/{sample}.2.trimmed.fastq.gz"
	resources:
		mem = 50
	threads: 1
	shell:
		"mv {input.f1} {output.f1} && mv {input.f2} {output.f2}" 

#Align reads to the reference genome using BWA
rule align_reads:
	input:
		"{dataset}/trimmed/{sample}.1.trimmed.fastq.gz",
		"{dataset}/trimmed/{sample}.2.trimmed.fastq.gz"
	output:
		"{dataset}/aligned/{sample}.bam"
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
		"{dataset}/aligned/{sample}.bam"
	output:
		"{dataset}/aligned/{sample}.sorted.bam"
	resources:
		mem = 8000
	threads: 4
	shell:
		"samtools sort -T {wildcards.dataset}/aligned/{wildcards.sample} -O bam -@ {threads} {input} > {output}"

#Index sorted bams
rule index_bams:
	input:
		"{dataset}/aligned/{sample}.sorted.bam"
	output:
		"{dataset}/aligned/{sample}.sorted.bam.bai"
	resources:
		mem = 50
	threads: 1
	shell:
		"samtools index {input}"

#Keep only properly paired reads from nuclear chromosomes (remove MT)
rule filter_properly_paired:
	input:
		bam = "{dataset}/aligned/{sample}.sorted.bam",
		index = "{dataset}/aligned/{sample}.sorted.bam.bai"
	output:
		"{dataset}/filtered/{sample}.filtered.bam"
	resources:
		mem = 100
	threads: 1
	params:
		chr_list = "1 10 11 12 13 14 15 16 17 18 19 2 20 21 22 3 4 5 6 7 8 9 X Y"
	shell:
		"samtools view -h -b -f 2 {input.bam} {params.chr_list} > {output}"

#Remove BWA entry from the BAM file header (conflicts with MarkDuplicates)
rule remove_bwa_header:
	input:
		"{dataset}/filtered/{sample}.filtered.bam"
	output:
		new_header = "{dataset}/filtered/{sample}.new_header.txt",
		bam = "{dataset}/filtered/{sample}.reheadered.bam"
	resources:
		mem = 100
	threads: 1
	shell:
		"samtools view -H {input} | grep -v 'ID:bwa' > {output.new_header} &&"
		"samtools reheader {output.new_header} {input} > {output.bam}"

#Remove duplicates using Picard MarkDuplicates
rule remove_duplicates:
	input:
		"{dataset}/filtered/{sample}.reheadered.bam"
	output:
		bam = protected("{dataset}/filtered/{sample}.no_duplicates.bam"),
		metrics = "{dataset}/metrics/{sample}.MarkDuplicates.txt"
	resources:
		mem = 2200
	threads: 4
	shell:
		"{config[picard_path]} MarkDuplicates I={input} O={output.bam} REMOVE_DUPLICATES=true METRICS_FILE= {output.metrics}"

#Count the number of reads per chromosome (QC metric)
rule reads_per_chromosome:
	input:
		"{dataset}/aligned/{sample}.sorted.bam"
	output:
		"{dataset}/metrics/{sample}.chr_counts.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"samtools view {input} | cut -f3 | sort | uniq -c > {output}"


#Sort BAM files by coordinates
rule sort_bams_by_name:
	input:
		"{dataset}/filtered/{sample}.no_duplicates.bam"
	output:
		"{dataset}/filtered/{sample}.sortedByName.bam"
	resources:
		mem = 8000
	threads: 4
	shell:
		"samtools sort -T {dataset}/filtered/{wildcards.sample} -O bam -@ {threads} -n {input} > {output}"

#Convert the bam file to a bed file of fragments
#First to bedpe and then to bed
rule bam_to_fragment_bed:
	input:
		"{dataset}/filtered/{sample}.sortedByName.bam"
	output:
		"{dataset}/bed/{sample}.fragments.bed.gz"
	resources:
		mem = 1000
	threads: 2
	shell:
		"bedtools bamtobed -bedpe -i {input} | python scripts/bedpe2bed.py --maxFragmentLength 1000 | sort -k 1,1 | gzip > {output}"

#Convert fragment bed files into cutsite bed files
rule fragments_to_cutsites_bed:
	input:
		"{dataset}/bed/{sample}.fragments.bed.gz"
	output:
		"{dataset}/bed/{sample}.cutsites.bed.gz"
	resources:
		mem = 1000
	threads: 2
	shell:
		"python scripts/fragmentsToCutSites.py --bed {input} | sort -k1,1 | gzip > {output}"

#Call peaks from cut site bed files
rule call_peaks:
	input:
		"{dataset}/bed/{sample}.cutsites.bed.gz"
	output:
		narrowPeak = "{dataset}/peaks/{sample}_peaks.narrowPeak"
	resources:
		mem = 1000
	threads: 1
	params:
		outdir = "{dataset}/peaks/",
		xls = "{dataset}/peaks/{sample}_peaks.xls"
	shell:
		"macs2 callpeak --nomodel -t {input} --shift 25 --extsize 50 -q 0.01 -n {wildcards.sample} --outdir {params.outdir} -f BED &&"
		"rm {params.xls}"

#Count the number of occurences of each fragment length in the fragments bed file
rule count_fragment_lengths:
	input:
		"{dataset}/bed/{sample}.fragments.bed.gz"
	output:
		"{dataset}/metrics/{sample}.fragment_lengths.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"zcat {input} | cut -f5 | sort -n | uniq -c > {output}"

#Convert fragment bed file into bigwig
rule convert_bed_to_bigwig:
	input:
		"{dataset}/bed/{sample}.fragments.bed.gz"
	output:
		bedgraph = "{dataset}/bigwig/{sample}.bg.gz",
		bigwig = "{dataset}/bigwig/{sample}.bw"
	params:
		bedgraph = "{dataset}/bigwig/{sample}.bg"
	resources:
		mem = 3000
	threads: 1
	shell:
		"bedtools genomecov -bga -i {input} -g {config[chromosome_lengths]} > {params.bedgraph} && "
		"bedGraphToBigWig {params.bedgraph} {config[chromosome_lengths]} {output.bigwig} && "
		"gzip {params.bedgraph}"

#Make sure that all final output files get created
rule make_all:
	input:
		expand("{dataset}/bigwig/{sample}.bw", sample=config["samples"]),
		expand("{dataset}/bed/{sample}.cutsites.bed.gz", sample=config["samples"]),
		expand("{dataset}/metrics/{sample}.fragment_lengths.txt", sample=config["samples"]),
		expand("{dataset}/metrics/{sample}.chr_counts.txt", sample=config["samples"])
	output:
		"{dataset}/out.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done' > {output}"

