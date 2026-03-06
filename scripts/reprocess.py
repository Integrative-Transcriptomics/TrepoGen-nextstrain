#!/usr/bin/env python3
"""Reprocess the contents of a Auspice JSON dataset file.
"""
import json, argparse
import pandas as pd

def parse_args():
	"""Parses command-line arguments for the dataset reprocessing script.

	Returns:
		argparse.Namespace: Parsed command-line arguments with the following attributes:
			dataset (str): Path to the input file (required).
			output (str): Path to the output file (required).
			metadata (str): Path to a metadata file (required).
			metadata-id-column (str): Column name in the metadata file to match IDs (required).
			metadata-add (list): List of metadata columns to add to the dataset.
			node-attr-remove (list): List of node attributes to remove from the dataset.
	"""
	parser = argparse.ArgumentParser(
		description="Reprocess the contents of a Auspice JSON dataset file.",
	)
	parser.add_argument("-d", "--dataset", type=str, required=True, help="The path to the input Auspice JSON dataset file.")
	parser.add_argument("-o", "--output", type=str, required=True, help="The path to the output Auspice JSON dataset file.")
	parser.add_argument("-m", "--metadata", type=str, required=True, help="The path to a metadata file to augment the dataset with.")
	parser.add_argument("--metadata-id-column", type=str, required=True, help="The column name in the metadata file to match IDs.")
	parser.add_argument("--metadata-add", type=str, nargs='*', help="The list of metadata columns to add to the dataset.")
	parser.add_argument("--node-attr-remove", type=str, nargs='*', help="The list of node attributes to remove from the dataset.")
	return parser.parse_args()

def main():
	""" Reprocess the contents of a Auspice JSON dataset file.

	Description:
		Reads an Auspice JSON dataset file, modifies its contents based on provided
		metadata and specified node attributes to remove, and writes the modified
		data back to the same file.

	Arguments:
		None. Arguments are parsed internally via parse_args().

	Returns:
		None.
	
	Side Effects:
		The input JSON file is modified in place, with specified node attributes removed
		and specified metadata columns added.
	"""
	# Parse command-line arguments.
	args = parse_args()

	# Load dataset.
	with open( args.dataset, 'r' ) as f:
		dataset = json.load( f )

	# Load metadata.
	meta_df = pd.read_csv( args.metadata, sep='\t' if args.metadata.endswith('.tsv') or args.metadata.endswith('.txt') else ',' )
	meta_df.set_index( args.metadata_id_column, inplace=True )

	# Reprocess dataset.
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
						if 'node_attrs' not in node :
							node[ 'node_attrs' ] = {}
						if meta_col not in node[ 'node_attrs' ] :
							node[ 'node_attrs' ][ meta_col ] = {}
						node[ 'node_attrs' ][ meta_col ][ 'value' ] = meta_df.at[ node_name, meta_col ]
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

	# Save dataset.
	with open( args.output, 'w' ) as f:
		json.dump( dataset, f )

if __name__ == '__main__':
	main()
