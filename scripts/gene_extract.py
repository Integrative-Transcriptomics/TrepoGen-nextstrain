#!/usr/bin/env python3
"""Extract gene sequences for a Nextstrain-TrepoGen gene dataset.

Description:
	Extracts gene-specific reference and sequences FASTA files
	using MUSIAL based on provided reference, annotation, and variant files.

Usage:
	python gene_extract.py
		-ir <reference.fasta>
		-ia <annotation.gff3>
		-ig <gene_name>
		-iv <variants.vcf>
		-or <output_reference.fasta>
		-os <output_sequences.fasta>
		-m <musial.jar>
		-q <root_sequence_identifier>
"""
import json
import subprocess
import argparse
import os
import tempfile
from pandas import DataFrame
from numpy import linspace
from textwrap import wrap
from glob import glob
from Bio import SeqIO
from Bio.Seq import Seq
from matplotlib import colormaps, colors

def parse_args():
	"""Parses command-line arguments for preparing data for a TrepoGen gene build using MUSIAL.

	Arguments:
		-ir, --reference (str, required): Reference genome sequence (FASTA) file.
		-ia, --annotation (str, required): Reference genome annotation (GFF3) file.
		-ig, --gene (str, required): The gene to extract sequence and annotation files for.
		-iv, --vcf (str, required): Variants (VCF) file.
		-or, --output-reference (str, required): Output path for the reference sequence file (FASTA).
		-os, --output-sequences (str, required): Output path for the sequences file (FASTA).
		-m, --musial (str, required): Path to the Musial JAR file or executable in PATH.
		-q, --root (str, required): Sequence identifier to use as the reference gene sequence.

	Returns:
		argparse.Namespace: Parsed command-line arguments with the following attributes:
		- reference (str): Path to the reference genome sequence file.
		- annotation (str): Path to the reference genome annotation file.
		- gene (str): Name of the gene to extract.
		- vcf (str): Path to the variants file.
		- output_reference (str): Path to the output reference sequence file.
		- output_sequences (str): Path to the output sequences file.
		- musial (str): Path to the Musial JAR file.
		- root (str): Sequence identifier to use as the reference gene sequence.
	"""
	parser = argparse.ArgumentParser(
		description="Prepare data for a TrepoGen gene build using MUSIAL.",
	)
	parser.add_argument("-ir", "--reference", type=str, required=True, help="Reference genome sequence (fasta) file.")
	parser.add_argument("-ia", "--annotation", type=str, required=True, help="Reference genome annotation (gff3) file.")
	parser.add_argument("-ig", "--gene", type=str, required=True, help="The gene to extract sequence and annotation files for.")
	parser.add_argument("-iv", "--vcf", type=str, required=True, help="Variants (vcf) file.")
	parser.add_argument("-or", "--output-reference", type=str, required=True, help="Output path for the reference sequence file (fasta).")
	parser.add_argument("-os", "--output-sequences", type=str, required=True, help="Output path for the sequences file (fasta).")
	parser.add_argument("-m", "--musial", type=str, required=True, help="Path to the Musial jar file.")
	parser.add_argument("-q", "--root", type=str, required=True, help="Sequence identifier to use as the reference gene sequence.")
	return parser.parse_args()

def write_config( reference : str, annotation : str, gene : str, vcf_input : str, tmpdir : str, min_coverage : int = 3, min_frequency : float = 0.8  ):
	"""Generates a MUSIAL configuration file for a Nextstrain-TrepoGen gene workflow.

	Description:
		Writes a features file and a configuration file for MUSIAL
		to be used in building a storage object for genomic sequences.

	Arguments:
		reference (str): Path to the reference genome file.
		annotation (str): Path to the genome annotation file.
		gene (str): Name of the gene to be processed.
		vcf_input (str): Path to the input VCF file.
		tmpdir (str): Directory where temporary files will be written.
		min_coverage (int, optional): Minimum coverage threshold. Defaults to 3.
		min_frequency (float, optional): Minimum frequency threshold. Defaults to 0.8.

	Returns:
		None.

	Side Effects:
		Writes a features file (features.tsv) and a configuration file (config.json) to the specified temporary directory.
	"""
	with open( f"{tmpdir}/features.tsv", "w+" ) as features_file :
		features_file.write( "id\tkey\tvalue\n" )
		features_file.write( f"{gene}\tgene\t{gene}\n" )
	config = {
		"minimalCoverage": min_coverage,
		"minimalFrequency": min_frequency,
		"maskFiltered": False,
		"skipAnnotation": True,
		"skipTyping": False,
		"reference": reference,
		"annotation": annotation,
		"features": f"{tmpdir}/features.tsv",
		"output": f"{tmpdir}/storage.json",
		"vcfFiles": [vcf_input]
	}
	with open( f"{tmpdir}/config.json", "w+" ) as config_file :
		json.dump( config, config_file )
		
def build_storage( musial_path : str, tmpdir : str ) :
	"""Builds the MUSIAL storage to extract genomic sequences from during a Nextstrain-TrepoGen gene workflow.

	Description:
		Runs MUSIAL build to create a storage object containing genomic sequences
		based on the provided configuration file.

	Arguments:
		musial_path (str): Path to the MUSIAL JAR file.
		tmpdir (str): Path to the temporary directory containing 'config.json'.

	Returns:
		None.

	Side Effects:
		Writes the storage object with information about genomic sequences to a file in the specified temporary directory.
		
	Raises:
		ChildProcessError: If the subprocess returns a non-zero exit code, indicating an error during execution.
	"""
	process = subprocess.Popen(
			[
				"java",
				"-Xmx4G",
				"-jar",
				musial_path,
				"build",
				"-C",
				f"{tmpdir}/config.json"
			],
			stdout=subprocess.PIPE,
			stderr=subprocess.PIPE,
		)
	stdout, stderr = process.communicate()
	stdout = stdout.decode(encoding="utf-8")
	stderr = stderr.decode(encoding="utf-8")
	if process.returncode != 0:
		raise ChildProcessError(f"Error during execution of subprocess.\nReturn code: {process.returncode}\n{stderr}")
	
def export_sequences( musial_path : str, gene : str, tmpdir : str ) :
	"""Exports gene sequences using MUSIAL during a Nextstrain-TrepoGen gene workflow.

	Description:
		Runs MUSIAL sequence to extract sequences for a specified gene feature
		and renames the output FASTA file to `sequences.fasta`.

	Arguments:
		musial_path (str): Path to the MUSIAL JAR file used for sequence extraction.
		gene (str): Name of the gene to extract sequences for.
		tmpdir (str): Path to a temporary directory for intermediate and output files.

	Returns:
		None.

	Side Effects:
		Writes a file named `sequences.fasta` in the specified temporary directory, containing the extracted gene sequences.

	Raises:
		ChildProcessError: If the Java subprocess reports an error (detected by "EXIT" in stderr).
		ValueError: If the number of `.fasta` files produced is not exactly one.
	"""
	process = subprocess.Popen(
			[
				"java",
				"-Xmx4G",
				"-jar",
				musial_path,
				"sequence",
				"--content",
				"nucleotide",
				"--locations",
				gene,
				"-I",
				f"{tmpdir}/storage.json",
				"--output",
				f"{tmpdir}"
			],
			stdout=subprocess.PIPE,
			stderr=subprocess.PIPE,
		)
	stdout, stderr = process.communicate()
	stdout = stdout.decode(encoding="utf-8")
	stderr = stderr.decode(encoding="utf-8")
	if process.returncode != 0:
		raise ChildProcessError( f"Error during execution of subprocess.\n{stderr}" )
	else :
		fasta_files = glob(f"{tmpdir}/*.fna")
		if len(fasta_files) != 1:
			raise ValueError(f"Expected exactly one .fna file, but found {len(fasta_files)}: {fasta_files}")
		os.rename(fasta_files[0], f"{tmpdir}/sequences.fasta")

def main():
	"""Extract gene specific data for a Nextstrain-TrepoGen gene dataset.

	Description:
	    Conducts the extraction of gene-specific reference and sequences FASTA files
		using MUSIAL based on provided reference, annotation, and variant files.

	Arguments:
		None. Arguments are parsed internally via `parse_args()`.
	
	Returns:
		None.

	Side Effects:
		Generates gene-specific reference and sequences FASTA files at the specified output paths.
	
	Raises:
		ValueError: If the reference sequence identifier cannot be parsed from the reference FASTA file.

	Note:
		This function assumes the existence of several helper functions and external dependencies,
		such as MUSIAL, Biopython's SeqIO, and pandas.
	"""
	
	args = parse_args()

	with tempfile.TemporaryDirectory( ) as tmpdir :
		# Write configuration files for MUSIAL.
		write_config( args.reference, args.annotation, args.gene, args.vcf, tmpdir )

		# Build storage for genomic sequences using the specified MUSIAL path and working directory.
		build_storage( args.musial, tmpdir )
	
		# Export sequences for the specified gene feature using MUSIAL.
		export_sequences( args.musial, args.gene, tmpdir )

		# Extract the seuence of the reference sample from the exported sequences and write to output file.
		samples_sequence_records = SeqIO.to_dict( SeqIO.parse( f"{tmpdir}/sequences.fasta", "fasta" ) )
		SeqIO.write( samples_sequence_records.values(), f"{args.output_sequences}", "fasta" )

		# Get the reference sequence record using the provided root identifier and write to output file.
		reference_seq_record = samples_sequence_records[ args.root ].seq
		with open( f"{args.output_reference}", "w+" ) as gene_reference_sequence_file :
			gene_reference_sequence_file.write( f">{args.root} [gene={args.gene}]\n" )
			gene_reference_sequence_file.write( "\n".join( wrap( str( reference_seq_record ), 60, break_on_hyphens=False ) ) + "\n" )

if __name__ == "__main__":
	main()
