#!/usr/bin/env python3
"""Reprocess the contents of a Auspice JSON dataset file.
"""
import json, argparse, yaml
import pandas as pd
from os.path import join, exists

# Define a set of red colors for resistance features.
REDS = [
	'#f86043',
	'#f5543c',
	'#f24834',
	'#ef3d2d',
	'#e63228',
	'#dd2924',
	'#d42120',
	'#cb181d',
	'#c1151b',
	'#b81319'
]

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

		If the source directory contains a resistance file (SOURCE/resistance.tsv), the script adds
		resistance information to the corresponding tree nodes.

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

	# Check if resistance information is available.
	resistance_file = join( "source", "data", args.source, 'resistance.tsv' )
	add_resistance_info = exists( resistance_file )
	if add_resistance_info :
		resistance_df = pd.read_csv( resistance_file, sep='\t' )
		feature_tags = { }
		def _get_resistance_info( node ) :
			# If node has children, skip - we only want to add resistance information to leaf nodes, i.e., actual samples.
			if node.get( 'children', [] ) :
				return
			
			# Set per resistance feature tags. Initially, all features are set to 'Sensitive'.
			node_tags = { }
			for feature in resistance_df[ "feature" ].unique().tolist() :
				node_tags[ feature ] = set( [ 'Sensitive' ] )

			# Get mutations per relevant gene.
			node_mutations = { }
			for gene in resistance_df[ "gene" ].unique().tolist() :
				node_mutations[gene] = set( node.get( 'branch_attrs', {} ).get( 'mutations', {} ).get( gene, [] ) )

			# Check if any resistance mutations are present for the node.
			for i, row in resistance_df.iterrows() :
				resistance_mutation = f"{row['ref']}{row['position']}{row['alt']}"
				if resistance_mutation in node_mutations.get( row['gene'], set() ) :
					node_tags[row['feature']].add( row['tag'] )

			# Construct final resistance labels for each feature based on tags.
			for feature, tags in node_tags.items() :
				if all( [ t == 'Sensitive' or t == 'Unknown' for t in tags ] ) :
					if 'Unknown' in tags :
						node_tags[ feature ] = 'Unknown'
					else :
						node_tags[ feature ] = 'Sensitive'
				else :
					node_tags[ feature ] = f"Resistant ({', '.join( sorted( [ tag for tag in tags if not tag in ['Unknown', 'Sensitive'] ] ) )})"
			# Update global feature tags with any new tags observed in the node.
			feature_tags.setdefault( feature, set() ).update( set( node_tags[ feature ] ) )
			return node_tags

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
		# Add resistance information if available.
		if add_resistance_info :
			resistance_info = _get_resistance_info( node )
			if resistance_info is not None :
				for feature, resistance_label in resistance_info.items() :
					if 'node_attrs' not in node :
						node[ 'node_attrs' ] = {}
					if feature not in node[ 'node_attrs' ] :
						node[ 'node_attrs' ][ feature ] = {}
					node[ 'node_attrs' ][ feature ][ 'value' ] = resistance_label
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

	# Add resistance feature colorings, if information was added.
	if add_resistance_info :
		for feature, tags in feature_tags.items() :
			coloring = {
				'key': feature,
				'title': " ".join( [ s.capitalize() for s in feature.split('_') ] ),
				'type': 'categorical',
			}
			scale = []
			# Set color for each tag. Resistant is set to reds, sensitive to blue, and unknown to grey.
			if 'Unknown' in tags :
				scale.append( [ 'Unknown', '#AAAAAA' ] )
			if 'Sensitive' in tags :
				scale.append( [ 'Sensitive', '#446DF6' ] )
			resistant_tags = [ t for t in tags if not t in ['Unknown', 'Sensitive'] ]
			if resistant_tags :
				global REDS
				for i, tag in enumerate( sorted( resistant_tags ) ) :
					color = REDS[ i % len(REDS) ]
					scale.append( [ tag, color ] )
			coloring['scale'] = scale
			colorings.append( coloring )

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
