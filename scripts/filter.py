#!/usr/bin/env python3
"""Remove sites from a VCF file based on a BED file.
"""
import argparse

def parse_args():
	"""Parses command-line arguments for the site removal script.

	Returns:
		argparse.Namespace: Parsed command-line arguments with the following attributes:
			input (str): Path to the input file (required).
			remove (str): Path to the BED file containing sites to remove (required).
			output (str): Path to the output file (required).
	"""
	parser = argparse.ArgumentParser(
		description="Remove sites specified in the input BED file from the input VCF file.",
	)
	parser.add_argument("-i", "--input", type=str, required=True, help="Path to the input auspice JSON dataset file.")
	parser.add_argument("-r", "--remove", type=str, required=True, help="Path BED file specifying sites to remove.")
	parser.add_argument("-o", "--output", type=str, required=True, help="Path to the output Auspice JSON dataset file.")
	return parser.parse_args()

def main():
	""" Remove sites specified in the input BED file from the input VCF file.

	Description:
		Modifies a VCF file by removing sites specified in a BED file. The output is a new VCF file with the specified sites removed.

	Arguments:
		None. Arguments are parsed internally via parse_args().

	Returns:
		None.
	"""
	# Parse command-line arguments.
	args = parse_args()

	# Load BED file.
	with open(args.remove, 'r') as f:
		bed_lines = f.readlines()
		bed_sites = set()
		for line in bed_lines:
			if line.startswith('#') or not line.strip():
				continue
			parts = line.strip().split('\t')
			if len(parts) < 3:
				continue
			chrom, start, end = parts[0], int(parts[1]), int(parts[2])
			for pos in range(start, end + 1):
				bed_sites.add((chrom, pos))

	# If no sites to remove, simply copy the input VCF to the output.
	if not bed_sites:
		with open(args.input, 'r') as f_in, open(args.output, 'w') as f_out:
			for line in f_in:
				f_out.write(line)
		return

	# Load VCF file.
	with open(args.input, 'r') as f:
		vcf_lines = f.readlines()

	# Filter VCF lines based on BED sites.
	with open(args.output, 'w') as f:
		for line in vcf_lines:
			if line.startswith('#'):
				f.write(line)
				continue
			parts = line.strip().split('\t')
			chrom, pos = parts[0], int(parts[1])
			if (chrom, pos) not in bed_sites:
				f.write(line)

if __name__ == '__main__':
	main()
