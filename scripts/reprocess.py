#!/usr/bin/env python3
"""Reprocess the contents of a Auspice JSON dataset file.

The script applies the following modifications to an Auspice JSON dataset file:
1. Adds CDS topology information, if files in SOURCE/topology/GENE.tsv are present for genes in the dataset.
2. Correctly tags date attributes as inferred, if the corresponding date in the metadata contains "
3. Adds gene coordinates from a genes.gff3 file, if it exists in the source directory.
"""
import json, argparse, yaml
import pandas as pd
from os.path import join, exists

def parse_args():
	"""Parses command-line arguments for the dataset reprocessing script.

	Returns:
		argparse.Namespace: Parsed command-line arguments with the following attributes:
			input (str): Path to the input file (required).
			output (str): Path to the output file (required).
			source (str): Path to the source directory (required). Expects to contain a meta.tsv file and optional topology tsv files in subdirectory topology/.
			id (str): (Optional) Column to use as node identifier for matching metadata to tree nodes. Defaults to "id".
	"""
	parser = argparse.ArgumentParser(
		description="Reprocess the contents of a Auspice JSON dataset file.",
	)
	parser.add_argument("-i", "--input", type=str, required=True, help="Path to the input auspice JSON dataset file.")
	parser.add_argument("-o", "--output", type=str, required=True, help="Path to the output Auspice JSON dataset file.")
	parser.add_argument("-s", "--source", type=str, required=True, help="Path to project specific source identifier.")
	parser.add_argument("--id", type=str, required=False, default="id", help="Column to use as node identifier for matching metadata to tree nodes.")
	return parser.parse_args()

def main():
	""" Reprocess the contents of a Auspice JSON dataset file.

	Description:
		Modifies an Auspice JSON dataset file by:
		- Adding CDS topology information, if files in SOURCE/topology/GENE.tsv are present for genes in the dataset.
		- Correctly tagging date attributes as inferred, if the corresponding date in the metadata contains "XX".

	Arguments:
		None. Arguments are parsed internally via parse_args().

	Returns:
		None.
	"""
	# Parse command-line arguments.
	args = parse_args()

	# Load dataset.
	with open( args.input, 'r' ) as f:
		dataset = json.load( f )

	# Load metadata.
	meta_file = join( "source", "data", args.source, 'meta.tsv' )
	# Validate that the file exists.
	if not exists( meta_file ) :
		raise FileNotFoundError( f"Metadata file '{meta_file}' not found." )
	meta_df = pd.read_csv( meta_file, sep='\t' if meta_file.endswith('.tsv') or meta_file.endswith('.txt') else ',' )
	if args.id not in meta_df.columns :
		raise ValueError( f"Specified ID column '{args.id}' not found in metadata file." )
	meta_df.set_index( args.id, inplace=True )

	# 1. Reprocess tree nodes to correctly tag date attributes as inferred.
	# Create a mapping of IDs to whether their date was inferred (contains "XX" or "/").
	date_inferred = { s: d for s, d in zip( meta_df.index.tolist(), meta_df.date.map( lambda d : "XX" in d or "/" in d ).tolist() ) }
	# Reprocess tree nodes.
	def _handle( node ) :
		node_name = node.get( 'name', None )
		# Adjust date inferred state.
		if 'num_date' in node.get( 'node_attrs', {} ) and node_name in date_inferred :
			node[ 'node_attrs' ][ 'num_date' ][ 'inferred' ] = date_inferred[ node_name ]
		# Recurse into children.
		for children in node.get( 'children', [] ) :
			_handle( children )
	# Start recursion at the tree root.
	_handle( dataset.get( 'tree', {} ) )

	# 2. If topology files are available, add topology information to matching genes and add description from source .yaml.
	feature_descriptions = {}
	with open ( join( "config", 'source.yaml' ), 'r' ) as f :
		source_config = yaml.safe_load( f )
		feature_descriptions = source_config.get( 'sources', {} ).get( args.source, {} ).get( 'genes', {} )
	for feature in dataset.get( 'meta', {} ).get( 'genome_annotations', {} ).keys() :
		if feature == 'nuc' :
			# Skip genome.
			continue
		topology_file = join( "source", "data", args.source, 'topology', f'{feature}.tsv' )
		if exists( topology_file ) :
			topology_df = pd.read_csv( topology_file, sep='\t', comment='#' )
			dataset.get( 'meta', {} ).get( 'genome_annotations', {} ).get( feature, {} )[ 'topology' ] = ",".join( [ f"{entry.type}:{entry.nuc_start}:{entry.nuc_end}" for entry in topology_df.itertuples() ] )
		if feature in feature_descriptions :
			dataset.get( 'meta', {} ).get( 'genome_annotations', {} ).get( feature, {} )[ 'info' ] = feature_descriptions[ feature ][ 'describe' ]

	# 3. Add gene coordinates from genes.gff3 file if it exists.
	# - Relevant genes are expected to have a 'gene_name' attribute in the GFF3 file.
	genes_file = join( "source", "data", args.source, 'genes.gff3' )
	# Check if the genes.gff3 file exists.
	if not exists( genes_file ) :
		print( f"Warning: Genes GFF3 file '{genes_file}' not found. Skipping `gene_coordinates` processing." )
	else :
		genes_df = pd.read_csv( genes_file, sep='\t', comment='#', header=None, names=['seqid','source','type','start','end','score','strand','phase','attributes'] )
		dataset[ 'meta' ].setdefault( 'extensions', {} )
		dataset[ 'meta' ][ 'extensions' ].setdefault( 'gene_coordinates', {} )
		for row_index, row in genes_df.iterrows() :
			attributes = { key: value for key, value in [ attr.split('=') for attr in row['attributes'].split(';') if '=' in attr ] }
			if 'gene_name' in attributes :
				gene_name = attributes['gene_name']
				start = int(row['start'])
				end = int(row['end'])
				strand = row['strand']
				dataset[ 'meta' ][ 'extensions' ][ 'gene_coordinates' ][ gene_name ] = { 'start': start, 'end': end, 'strand': strand }
			
	# Save the reprocessed dataset.
	with open( args.output, 'w' ) as f:
		json.dump( dataset, f )

if __name__ == '__main__':
	main()
