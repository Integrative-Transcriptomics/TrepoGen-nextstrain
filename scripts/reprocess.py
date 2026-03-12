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
			metadata-id-column (str): Column name in the metadata file to match IDs (required).
			metadata-add (list): List of metadata columns to add to the dataset.
			node-attr-remove (list): List of node attributes to remove from the dataset.
	"""
	parser = argparse.ArgumentParser(
		description="Reprocess the contents of a Auspice JSON dataset file.",
	)
	parser.add_argument("-i", "--input", type=str, required=True, help="Path to the input auspice JSON dataset file.")
	parser.add_argument("-o", "--output", type=str, required=True, help="Path to the output Auspice JSON dataset file.")
	parser.add_argument("-s", "--source", type=str, required=True, help="Path to project specific source identifier.")
	parser.add_argument("--metadata-id-column", type=str, required=True, help="Column name in the metadata file to match IDs.")
	parser.add_argument("--metadata-add", type=str, nargs='*', help="List of metadata columns to add to the dataset.")
	parser.add_argument("--node-attr-remove", type=str, nargs='*', help="List of node attributes to remove from the dataset.")
	return parser.parse_args()

def main():
	""" Reprocess the contents of a Auspice JSON dataset file.

	Description:
		Reads an Auspice JSON dataset file, modifies its contents based on provided
		metadata and specified node attributes to remove, and writes the modified
		data back to the specified output file.

		If the source directory contains topology files for genes (SOURCE/topology/GENE.tsv),
		the script adds topology information to the corresponding gene annotations in the dataset.
		Also adds additional descriptions for features if available from the source.yaml configuration.

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
	meta_df.set_index( args.metadata_id_column, inplace=True )

	# Reprocess tree nodes.
	def _handle( node ) :
		# Remove specified node attributes.
		for attr in ( args.node_attr_remove if args.node_attr_remove is not None else [] ) :
			del node.get( 'node_attrs', {} )[ attr ]

		# Add specified metadata columns.
		if 'node_attrs' in node :
			node_name = node.get( 'name', None )
			if node_name is not None and node_name in meta_df.index :
				for meta_col in ( args.metadata_add if args.metadata_add is not None else [] ) :
					if meta_col in meta_df.columns :
						meta_label = str(meta_col).capitalize()
						if 'node_attrs' not in node :
							node[ 'node_attrs' ] = {}
						if meta_label not in node[ 'node_attrs' ] :
							node[ 'node_attrs' ][ meta_label ] = {}
						node[ 'node_attrs' ][ meta_label ][ 'value' ] = meta_df.at[ node_name, meta_col ]

		# Recurse into children.
		for children in node.get( 'children', [] ) :
			_handle( children )

	# Start recursion at the tree root.
	_handle( dataset.get( 'tree', {} ) )

	# Remove specified colorings of attributes that were removed.
	colorings = dataset.get( 'meta', {} ).get( 'colorings', [] )
	for attr in ( args.node_attr_remove if args.node_attr_remove is not None else [] ) :
		for coloring in colorings :
			if coloring.get( 'key', None ) == attr :
				colorings.remove( coloring )

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
