#Read barcode sequecnes from disk
barcode1 = open(snakemake.input["barcode1"]).readline().rstrip()
barcode2 = open(snakemake.input["barcode1"]).readline().rstrip()

#Construct the full i7 and i5 primer sequences
primer1 = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC" + barcode1 + "ATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAA"
primer2 = "CTGTCTCTTATACACATCTGACGCTGCCGACGA" + barcode2 + "GTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAA"

#Construct skewer command
output_prefix = snakemake.params["prefix"]
skewer_command = " ".join(["skewer -x", primer1, "-y", primer2, "-m pe", 
	snakemake.input["fastq1"], snakemake.input["fastq2"], "-z -o", output_prefix])
print(skewer_command)