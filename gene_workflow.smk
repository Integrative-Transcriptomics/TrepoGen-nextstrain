configfile: "config/global.yaml"
configfile: "config/source.yaml"
configfile: "config/gene.yaml"

# genes: List of genes to use in the workflow (cf. config/genes.yaml); a dataset will be generated for each.
genes = list(config.get("genes", {}).keys())

# sources: List of source datasets to use in the workflow; expected to match (VAL) in source/data/<VAL>.
sources = ["TPASS-308"]

# subsets: List of variant subsets to use in the workflow; expected to match (VAL) in source/data/*/variants/<VAL>.vcf.
subsets = ["snv-indel"]

# source_files: List of files that are part of the source data; files are expected to be present in the source/data/{source} directory.
source_files = [
	"sequence.fasta", # Reference genome sequence in FASTA format.
	"annotation.gff3", # Annotation file in GFF3 format.
	"meta.csv", # Metadata file in CSV format containing information about the single samples.
	"meta_colors.tsv", # TSV file defining colors for the metadata.
	"meta_configuration.json", # Auspice configuration file in JSON format providing settings for the meta data of the dataset.
]

# work_files: List of files that are generated during the workflow; files will be generated in the .work/{source}_{subset}_{gene}/ directory.
work_files = [
	"reference.fasta", # Reference gene sequence in FASTA format.
	"annotation.gff3", # Annotation file in GFF3 format for the gene.
	"sequences.fasta", # Sample sequences in FASTA format for the gene.
	"gene_features.txt", # Text file containing names of features in the gene.
	"gene_meta.csv", # Types file in CSV format containing information about selected feature types, e.g. extracellular domains.
	"gene_configuration.json", # Auspice configuration file in JSON format for the selected feature types.
	"meta.csv", # Metadata file in CSV format containing information about the samples.
	"initial.nwk", # Initial tree in Newick format (NWK) before refinement.
	"tree.nwk", # Refined tree in Newick format (NWK) after refinement.
	"traits.json", # Traits in JSON format for the tree.
	"branch_lengths.json", # Branch lengths in JSON format for the tree.
	"nucleotide_mutations.json", # Nucleotide mutations in JSON format for the tree.
	"amino_acid_mutations.json", # Amino acid mutations in JSON format for the tree.
	"description.md", # Description file in Markdown format (MD) for the dataset.
]

# Defines the directory structure and I/O files for each source, subset, and gene.
rule all:
	input:
		expand("source/data/{source}/{source_file}", source=sources, source_file=source_files),
		expand("source/data/{source}/variants/{subset}.vcf", source=sources, subset=subsets),
		expand("source/data/{source}/variants/{subset}.tsv", source=sources, subset=subsets),
		expand(".work/{source}_{subset}_{gene}/filtered.vcf", source=sources, subset=subsets, gene=genes),
		expand(".work/{source}_{subset}_{gene}/{work_file}", source=sources, subset=subsets, gene=genes, work_file=work_files),
		expand("datasets/{source}_{subset}_{gene}.json", source=sources, subset=subsets, gene=genes),
		expand("datasets/{source}_{subset}_{gene}_reprocessed.json", source=sources, subset=subsets, gene=genes),
		"source/misc/gene_display_configuration.json", # Independent file for gene dataset display defaults.
		"source/geo/color.tsv", # Independent file for geo colors.
		"source/geo/loc.tsv" # Independent file for geo coordinates.

# Generates an index of the input variants.
rule index:
	input:
		"source/data/{source}/variants/{subset}.vcf",
	output:
		"source/data/{source}/variants/{subset}.tsv",
	shell:
		"""
		augur index \
		  --sequences {input} \
		  --output {output}
		"""

# Filters variants (optional) based on metadata and excludes specified samples. If no filter is defined, it simply passes the data through.
rule filter:
	input:
		variants="source/data/{source}/variants/{subset}.vcf",
		index="source/data/{source}/variants/{subset}.tsv",
		metadata="source/data/{source}/meta.csv",
	output:
		".work/{source}_{subset}_{gene}/filtered.vcf",
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

# Prepares a gene build by processing the reference genome, annotation, and variants into gene specific data.
rule extract:
	input:
		reference="source/data/{source}/sequence.fasta",
		annotation="source/data/{source}/annotation.gff3",
		variants=rules.filter.output,
	output:
		reference=".work/{source}_{subset}_{gene}/reference.fasta",
		annotation=".work/{source}_{subset}_{gene}/annotation.gff3",
		features=".work/{source}_{subset}_{gene}/gene_features.txt",
		sequences=".work/{source}_{subset}_{gene}/sequences.fasta",
		gene_meta=".work/{source}_{subset}_{gene}/gene_meta.csv",
		gene_configuration=".work/{source}_{subset}_{gene}/gene_configuration.json",
	params:
		musial=config.get("musial", "musial"),
		root=lambda wc: config["sources"][wc.source].get("reference_sample", ""),
	shell:
		"""
		python scripts/gene_extract.py \
			-ir {input.reference} \
			-ia {input.annotation} \
			-iv {input.variants} \
			-ig {wildcards.gene} \
			-or {output.reference} \
			-oa {output.annotation} \
			-og {output.features} \
			-os {output.sequences} \
			-ot {output.gene_meta} \
			-oc {output.gene_configuration} \
			-m {params.musial} \
			-q {params.root}
		"""

# Merges metadata files of the source data and the gene build.
rule metadata:
	input:
		metadata_source="source/data/{source}/meta.csv",
		metadata_gene=rules.extract.output.gene_meta,
	output:
		".work/{source}_{subset}_{gene}/meta.csv"
	params:
		metadata_id=lambda wc: config["sources"][wc.source].get("meta_identifier", "name strain id"),
	shell:
		"""
		augur merge \
			--metadata msource={input.metadata_source} mgene={input.metadata_gene} \
			--metadata-id-columns {params.metadata_id} \
			--output-metadata {output}
		"""

# Concatenates color definition files of the source data as needed for the dataset.
rule colors:
	input:
		meta_colors="source/data/{source}/meta_colors.tsv",
		geo_colors="source/geo/color.tsv",
	output:
		".work/{source}_{subset}_{gene}/colors.tsv",
	shell:
		"""
		cat {input.meta_colors} {input.geo_colors} >> {output}
		"""

# Builds the phylogenetic tree from the sequence alignment.
rule tree:
	input:
		alignment=rules.extract.output.sequences,
	output:
		".work/{source}_{subset}_{gene}/initial.nwk",
	params:
		method_cl=lambda wc: config["genes"][wc.gene].get("tree", {}).get("method_cl", "--method iqtree"),
	shell:
		"""
		augur tree --alignment {input.alignment} \
			--output {output} \
			{params.method_cl}
		"""

# Refines the phylogenetic tree; performs date inference and branch length estimation.
rule refine:
	input:
		alignment=rules.extract.output.sequences,
		metadata=rules.metadata.output,
		tree=rules.tree.output,
	output:
		tree=".work/{source}_{subset}_{gene}/tree.nwk",
		branch_lengths=".work/{source}_{subset}_{gene}/branch_lengths.json",
	params:
		seed=config.get("seed", 1),
		iterations=config.get("refine", {}).get("iterations", 1),
		precision=config.get("refine", {}).get("precision", 1),
		metadata_id=lambda wc: config["sources"][wc.source].get("meta_identifier", "name strain id"),
		clock_rate_cl=lambda wc: config["genes"][wc.gene].get("refine", {}).get("clock_rate_cl", ""),
		year_bounds_cl=lambda wc: config["sources"][wc.source].get("refine", {}).get("year_bounds_cl", ""),
	shell:
		"""
		augur refine --tree {input.tree} \
			--alignment {input.alignment} \
			--metadata {input.metadata} \
			--metadata-id-columns {params.metadata_id} \
			--timetree \
			--use-fft \
			--max-iter {params.iterations} \
			--covariance \
			--keep-polytomies \
			--precision {params.precision} \
			--coalescent opt \
			--date-format %Y-%m-%d \
			--date-inference marginal \
			--date-confidence \
			--branch-length-inference marginal \
			--seed {params.seed} \
			--output-tree {output.tree} \
			--output-node-data {output.branch_lengths} \
			{params.clock_rate_cl} \
			{params.year_bounds_cl}
		"""

# Computes the ancestral sequences for the phylogenetic tree, inferring nucleotide mutations.
rule ancestral:
	input:
		alignment=rules.extract.output.sequences,
		reference=rules.extract.output.reference,
		tree=rules.refine.output.tree,
	output:
		".work/{source}_{subset}_{gene}/nucleotide_mutations.json",
	params:
		seed=config.get("seed", 1),
	shell:
		"""
		augur ancestral --tree {input.tree} \
			--alignment {input.alignment} \
			--root-sequence {input.reference} \
			--inference marginal \
			--keep-ambiguous \
			--keep-overhangs \
			--seed {params.seed} \
			--output-node-data {output}
		"""

# Computes the ancestral traits for the phylogenetic tree, inferring traits from the metadata.
rule traits:
	input:
		tree=rules.refine.output.tree,
		metadata=rules.metadata.output,
	output:
		".work/{source}_{subset}_{gene}/traits.json",
	params:
		metadata_id=lambda wc: config["sources"][wc.source].get("meta_identifier", "name strain id"),
		columns=lambda wc: config["sources"][wc.source].get("traits", {}).get("columns", "country date"),
	shell:
		"""
		augur traits --tree {input.tree} \
			--metadata {input.metadata} \
			--metadata-id-columns {params.metadata_id} \
			--columns {params.columns} \
			--output-node-data {output}
		"""

# Translates nucleotide sequences of specified gene features into amino acid sequences.
rule translate:
	input:
		tree=rules.refine.output.tree,
		sequences=rules.ancestral.output,
		annotation=rules.extract.output.annotation,
	output:
		".work/{source}_{subset}_{gene}/amino_acid_mutations.json",
	params:
		genes=rules.extract.output.features,
	run:
		shell(
			"""
			augur translate \
				--tree {input.tree} \
				--ancestral-sequences {input.sequences} \
				--reference-sequence {input.annotation} \
				--output-node-data {output} \
				--genes {params.genes}
			"""
		)
		shell(
			"""
			python scripts/gene_color_annotations.py -f {output}
			"""
		)

# Generates a description file for the dataset.
rule describe:
	output:
		".work/{source}_{subset}_{gene}/description.md",
	params:
		content=lambda wc: "\n".join([
			config["describe"]["gene"][wc.gene],
			config["describe"]["gene"]["postscript"],
			config["describe"]["source"][wc.source],
			config["describe"]["subset"][wc.subset],
			"---",
			config["describe"]["resources"],
			config["describe"]["background"],
			config["describe"]["funding"]
		]),
	shell:
		"""
		echo '{params.content}' > {output}
		"""

# Exports the refined tree and associated data into a format suitable for visualization in Auspice.
rule export:
	input:
		tree=rules.refine.output.tree,
		metadata=rules.metadata.output,
		colors=rules.colors.output,
		meta_config="source/data/{source}/meta_configuration.json",
		gene_config=rules.extract.output.gene_configuration,
		display_config="source/misc/gene_display_configuration.json",
		description=rules.describe.output,
		coordinates="source/geo/loc.tsv",
		branch_lengths=rules.refine.output.branch_lengths,
		traits=rules.traits.output,
		nucleotide_mutations=rules.ancestral.output,
		amino_acid_mutations=rules.translate.output,
	output:
		"datasets/{source}_{subset}_{gene}.json",
	params:
		metadata_id=lambda wc: config["sources"][wc.source].get("meta_identifier", "name strain id"),
		title=lambda wc: f"'{config.get('export', {}).get('title', 'undefined')} ({wc.gene})'",
		maintainers=config.get("export", {}).get("maintainers", "undefined"),
		build_url=config.get("export", {}).get("build_url", "undefined"),
	shell:
		"""
		augur export v2 \
			--tree {input.tree} \
			--metadata {input.metadata} \
			--metadata-id-columns {params.metadata_id} \
			--node-data {input.branch_lengths} {input.traits} {input.nucleotide_mutations} {input.amino_acid_mutations} \
			--auspice-config {input.meta_config} {input.gene_config} {input.display_config} \
			--title {params.title} \
			--maintainers {params.maintainers} \
			--build-url {params.build_url} \
			--description {input.description} \
			--colors {input.colors} \
			--lat-longs {input.coordinates} \
			--output {output} \
			--include-root-sequence-inline
		"""

# Reprocesses the exported dataset to add or remove specific metadata and node attributes.
rule reprocess:
	input:
		dataset=rules.export.output,
		metadata="source/data/{source}/meta.csv",
	output:
		"datasets/{source}_{subset}_{gene}_reprocessed.json",
	params:
		metadata_id=lambda wc: config["sources"][wc.source].get("meta_identifier", "name strain id"),
		metadata_add="--metadata-add label",
		node_attr_remove="--node-attr-remove resistance_mutations",
	shell:
		"""
		python scripts/reprocess_dataset.py \
			--dataset {input.dataset} \
			--output {output} \
			--metadata {input.metadata} \
			--metadata-id-column {params.metadata_id} \
			{params.metadata_add} \
			{params.node_attr_remove}
		"""
