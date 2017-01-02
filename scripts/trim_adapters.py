#Read barcode sequecnes from disk
barcode1 = open(snakemake.input["barcode1"]).readline().rstrip()
barcode2 = open(snakemake.input["barcode1"]).readline().rstrip()
print(barcode1)
print(barcode2)