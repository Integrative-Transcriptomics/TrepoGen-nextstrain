# Load configuration files for the workflow.
configfile: "config/global.yaml" # Global configurations.
configfile: "config/source.yaml" # Data-source specific configurations.

sources = ["TPA-SS14"]

# List of source files that are part of the source data; files are expected to be present in the source/data/{source} directory.
source_files = [
	"sequence.fasta", # Reference genome sequence in FASTA format.
	"annotation.gff3", # Annotation file in GFF3 format.
	"variants.vcf", # VCF file containing the variants for the samples in the dataset.
	"mask.bed", # (Optional empty) BED file for masking regions in the genome for tree building.
	"meta.csv", # Metadata file in CSV format containing information about the single samples.
	"meta_colors.tsv", # TSV file defining colors for the metadata.
	"meta_configuration.json", # Auspice configuration file in JSON format providing settings for the metadata of the dataset.
]

# ! The list of genes to use have to be specified in config/source.yaml under the {source} > genes section.

# List of work files that are generated during the workflow; files will be generated in the .work/{source}/ directory.
work_files = [
	"initial.nwk", # Initial tree in Newick format (NWK) before refinement.
	"tree.nwk", # Refined tree in Newick format (NWK) after refinement.
	"branch_lengths.json", # Branch lengths in JSON format for the tree.
	"nucleotide_mutations.json", # Nucleotide mutations in JSON format for the tree.
	"traits.json", # Traits in JSON format for the tree.
	"amino_acid_mutations.json", # Amino acid mutations in JSON format for the tree.
	"description.md", # Description file in Markdown format (MD) for the dataset.
]

# Defines the directory structure and I/O files for each source, subset, and gene.
rule all:
	input:
		expand("source/data/{source}/{source_file}", source=sources, source_file=source_files),
		expand("source/data/{source}/variants.tsv", source=sources), # Sample name index file created by augur.
		expand(".work/{source}/filtered.vcf", source=sources), # Filtered variants after applying metadata-based filters.
		expand(".work/{source}/{work_file}", source=sources, work_file=work_files),
		expand("datasets/{source}_raw.json", source=sources), # Raw dataset exported by augur export before reprocessing.
		expand("datasets/{source}.json", source=sources), # Final dataset after reprocessing.
		"source/misc/display_configuration.json", # Source independent file for dataset display defaults.
		"source/geo/color.tsv", # Source independent file for geo colors.
		"source/geo/loc.tsv" # Source independent file for geo coordinates.

# Generates an index of the input variants.
rule index:
	input:
		"source/data/{source}/variants.vcf",
	output:
		"source/data/{source}/variants.tsv", # Index is stored in source directory for reuse.
	shell:
		"""
		augur index \
		  --sequences {input} \
		  --output {output}
		"""

# Filters variants (optional) based on metadata and excludes specified samples. If no filter is defined, it simply passes the data through.
rule filter:
	input:
		variants="source/data/{source}/variants.vcf",
		index="source/data/{source}/variants.tsv",
		metadata="source/data/{source}/meta.csv",
	output:
		".work/{source}/filtered.vcf",
	params:
		metadata_id=lambda wc: config["sources"][wc.source].get("meta_identifier", "name strain id"),
		query_cl=lambda wc: config["sources"][wc.source].get("filter", {}).get("query_cl", ""),  # Optional query command line argument in config to filter variants.
	run:
		if bool(params.query_cl):
			shell(
				"""
				augur filter --sequences {input.variants} \
					--sequence-index {input.index} \
					--metadata {input.metadata} \
					--metadata-id-columns {params.metadata_id} \
					--output-sequences {output} \
					{params.query_cl}
				"""
			)
		else:
			# If no filter is defined, just pass the file through.
			shell("cp {input.variants} {output}")

# Concatenates color definition files of the source data as needed for the dataset.
rule colors:
	input:
		meta_colors="source/data/{source}/meta_colors.tsv",
		geo_colors="source/geo/color.tsv",
	output:
		".work/{source}/colors.tsv",
	shell:
		"""
		cat {input.meta_colors} {input.geo_colors} >> {output}
		"""

# Builds the initial phylogenetic tree from the variants.
rule tree:
	input:
		variants=rules.filter.output,
		reference="source/data/{source}/sequence.fasta",
		mask="source/data/{source}/mask.bed",
	output:
		".work/{source}/initial.nwk",
	params:
		method_cl=lambda wc: config["sources"][wc.source].get("tree", {}).get("method_cl", "--method iqtree"),
	shell:
		"""
		augur tree --alignment {input.variants} \
			--vcf-reference {input.reference} \
			--exclude-sites {input.mask} \
			--output {output} \
			{params.method_cl}
		"""

# Refines the phylogenetic tree; performs date inference and branch length estimation.
rule refine:
	input:
		variants=rules.filter.output,
		reference="source/data/{source}/sequence.fasta",
		metadata="source/data/{source}/meta.csv",
		tree=rules.tree.output,
	output:
		tree=".work/{source}/tree.nwk",
		branch_lengths=".work/{source}/branch_lengths.json",
	params:
		seed=config.get("seed", 1),
		iterations=lambda wc: config["sources"][wc.source].get("refine", {}).get("iterations", 1),
		precision=lambda wc: config["sources"][wc.source].get("refine", {}).get("precision", 1),
		metadata_id=lambda wc: config["sources"][wc.source].get("meta_identifier", "name strain id"),
		clock_rate_cl=lambda wc: config["sources"][wc.source].get("refine", {}).get("clock_rate_cl", ""),
		root_cl=lambda wc: config["sources"][wc.source].get("refine", {}).get("root_cl", ""),
		year_bounds_cl=lambda wc: config["sources"][wc.source].get("refine", {}).get("year_bounds_cl", ""),
	shell:
		"""
		augur refine --tree {input.tree} \
			--alignment {input.variants} \
			--vcf-reference {input.reference} \
			--metadata {input.metadata} \
			--metadata-id-columns {params.metadata_id} \
			--timetree \
			--use-fft \
			--max-iter {params.iterations} \
			--covariance \
			--keep-polytomies \
			--precision {params.precision} \
			--coalescent 'opt' \
			--date-format '%Y-%m-%d' \
			--date-inference 'joint' \
			--date-confidence \
			--branch-length-inference 'joint' \
			--seed {params.seed} \
			--output-tree {output.tree} \
			--output-node-data {output.branch_lengths} \
			--verbosity 4 \
			{params.clock_rate_cl} \
			{params.root_cl} \
			{params.year_bounds_cl}
		"""

# Computes the ancestral sequences for the phylogenetic tree, inferring nucleotide mutations.
rule ancestral:
	input:
		variants=rules.filter.output,
		reference="source/data/{source}/sequence.fasta",
		tree=rules.refine.output.tree,
	output:
		node_data=".work/{source}/nucleotide_mutations.json",
		ancestral_variants=".work/{source}/nucleotide_mutations.vcf",
	params:
		seed=config.get("seed", 1),
	shell:
		"""
		augur ancestral --tree {input.tree} \
			--alignment {input.variants} \
			--vcf-reference {input.reference} \
			--inference marginal \
			--keep-ambiguous \
			--keep-overhangs \
			--seed {params.seed} \
			--output-node-data {output.node_data} \
			--output-vcf {output.ancestral_variants}
		"""

# Computes the ancestral traits for the phylogenetic tree, inferring traits from the metadata.
rule traits:
	input:
		tree=rules.refine.output.tree,
		metadata="source/data/{source}/meta.csv",
	output:
		".work/{source}/traits.json",
	params:
		metadata_id=lambda wc: config["sources"][wc.source].get("meta_identifier", "name strain id"),
		columns=lambda wc: config["sources"][wc.source].get("traits", {}).get("columns", "date"),
	shell:
		"""
		augur traits --tree {input.tree} \
			--metadata {input.metadata} \
			--metadata-id-columns {params.metadata_id} \
			--columns {params.columns} \
			--output-node-data {output}
		"""

# Translates nucleotide sequences of specified genes into amino acid sequences.
rule translate:
	input:
		tree=rules.refine.output.tree,
		ancestral_variants=rules.ancestral.output.ancestral_variants,
		reference="source/data/{source}/sequence.fasta",
		annotation="source/data/{source}/annotation.gff3",
	output:
		".work/{source}/amino_acid_mutations.json",
	params:
		genes=lambda wc: " ".join(config["sources"][wc.source].get("genes", {}).keys()).strip(), # List of genes to translate, expected to match the gene names in the annotation.gff3 file.
	run:
		if len(params.genes) > 0:
			shell(
				"""
				augur translate \
					--tree {input.tree} \
					--ancestral-sequences {input.ancestral_variants} \
					--reference-sequence {input.annotation} \
					--vcf-reference {input.reference} \
					--output-node-data {output} \
					--genes {params.genes}
				"""
			)
		else:
			# If no genes are defined, create an empty dummy file.
			shell(
				"""
				echo '{{}}' > {output}
				"""
			)

# Generates a description file for the dataset.
rule describe:
	output:
		".work/{source}/description.md",
	params:
		content=lambda wc: "\n".join([
			config["describe"]["preface"],
			config["sources"][wc.source]["describe"],
			"---",
			config["describe"]["postscript"]
		]),
	shell:
		"""
		echo "{params.content}" > {output}
		"""

# Exports the refined tree and associated data into a format suitable for visualization in Auspice.
rule export:
	input:
		tree=rules.refine.output.tree,
		metadata="source/data/{source}/meta.csv",
		colors=rules.colors.output,
		meta_config="source/data/{source}/meta_configuration.json",
		display_defaults="source/misc/display_configuration.json",
		description=rules.describe.output,
		coordinates="source/geo/loc.tsv",
		branch_lengths=rules.refine.output.branch_lengths,
		traits=rules.traits.output,
		nucleotide_mutations=rules.ancestral.output.node_data,
		amino_acid_mutations=rules.translate.output,
	output:
		"datasets/{source}_raw.json",
	params:
		metadata_id=lambda wc: config["sources"][wc.source].get("meta_identifier", "name strain id"),
		metadata_columns=lambda wc: config["sources"][wc.source].get("export", {}).get("metadata_columns", ""),
		title=lambda wc: f"'{config.get('export', {}).get('title', 'Nextstrain Dataset')}'",
		maintainers=config.get("export", {}).get("maintainers", "None"),
		build_url=config.get("export", {}).get("build_url", "Unknown"),
	shell:
		"""
		augur export v2 \
			--tree {input.tree} \
			--metadata {input.metadata} \
			--metadata-columns {params.metadata_columns} \
			--metadata-id-columns {params.metadata_id} \
			--node-data {input.branch_lengths} {input.traits} {input.nucleotide_mutations} {input.amino_acid_mutations} \
			--auspice-config {input.meta_config} {input.display_defaults} \
			--title {params.title} \
			--maintainers {params.maintainers} \
			--build-url {params.build_url} \
			--description {input.description} \
			--colors {input.colors} \
			--lat-longs {input.coordinates} \
			--output {output} \
			--include-root-sequence-inline
		"""

# Reprocesses the exported dataset.
rule reprocess:
	input:
		dataset=rules.export.output,
	output:
		"datasets/{source}.json",
	params:
		source=lambda wc: wc.source,
		metadata_id=lambda wc: config["sources"][wc.source].get("meta_identifier", "strain"),
	shell:
		"""
		python scripts/reprocess.py \
			--input {input.dataset} \
			--output {output} \
			--source {params.source} \
			--id {params.metadata_id}
		"""
