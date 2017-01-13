import os
import sys
import argparse
import fileinput
import subprocess
import gzip

parser = argparse.ArgumentParser(description = "Extract Nextera cut sites from fragments bed file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--bed", help = "Path to the fragments bed file.")
args = parser.parse_args()

bed_file = gzip.open(args.bed)
for line in bed_file:
	line = line.rstrip()
	fields = line.split("\t")
	new_start = int(fields[1]) + 4
	new_end = int(fields[2]) -4
	new_length = int(fields[4]) - 7
	cut1 = "\t".join([fields[0], str(new_start), str(new_start), fields[3], str(new_length), "+"])
	cut2 = "\t".join([fields[0], str(new_end), str(new_end), fields[3], str(new_length), "+"])
	print(cut1)
	print(cut2)