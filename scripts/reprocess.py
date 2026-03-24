#!/usr/bin/env python3
"""Reprocess the contents of a Auspice JSON dataset file.
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
			source (str): Path to the source directory (required). Expects to contain a meta.csv file and optional topology tsv files in subdirectory topology/.
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
	meta_file = join( "source", "data", args.source, 'meta.csv' )
	meta_df = pd.read_csv( meta_file, sep='\t' if meta_file.endswith('.tsv') or meta_file.endswith('.txt') else ',' )
	if args.id not in meta_df.columns :
		raise ValueError( f"Specified ID column '{args.id}' not found in metadata file." )
	meta_df.set_index( args.id, inplace=True )

	# Create a mapping of IDs to whether their date is inferred (contains "XX").
	date_inferred = { s: d for s, d in zip( meta_df.index.tolist(), meta_df.date.map( lambda d : "XX" in d ).tolist() ) }

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

	# If topology files are present, add topology information to matching genes; also add additional description if available from source .yaml.
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

	# Save dataset.
	with open( args.output, 'w' ) as f:
		json.dump( dataset, f )

if __name__ == '__main__':
	main()
