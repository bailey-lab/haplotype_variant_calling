configfile: 'call_variants.yaml'

rule all:
	input:
		#seekdeep_tables=config['output_folder']+'/allSelectedClustersInfo.tab.txt.gz'
		#called_aa_variants=config['output_folder']+'/variantCallOnSeqAndProtein/reports/AAChangesInfo.tsv.gz'
		prevalence_table=config['output_folder']+'/amino_acid_prevalences.tsv'

rule concatenate_seekdeep_tables:
	'''
	concatenates seekdeep tables into a single table suitable for variant calling
	'''
	input:
		seekdeep_folder=config['seekdeep_popclustering_folder']
	output:
		seekdeep_tables=config['output_folder']+'/allSelectedClustersInfo.tab.txt.gz'
	script:
		'scripts/concatenate_seekdeep_tables.py'

rule call_variants:
	'''
	runs Nick Hathaway's variantCallOnSeqAndProtein script to call variants
	'''
	input:
		seekdeep_tables=config['output_folder']+'/allSelectedClustersInfo.tab.txt.gz',
		genome=config['genome_file'],
		gff=config['gff_file'],
		targets_of_interest=config['seekdeep_targets_file']
	params:
		bindings=config['bindings'],
		sif_file=config['sif_file'],
		output_folder=config['output_folder']+'/variantCallOnSeqAndProtein'
	threads: 14
	output:
		called_aa_variants=config['output_folder']+'/variantCallOnSeqAndProtein/reports/AAChangesInfo.tsv.gz'
	shell:
		'''
		singularity exec {params.bindings} {params.sif_file} SeekDeep \
		variantCallOnSeqAndProtein --resultsFnp {input.seekdeep_tables} \
		--genomeFnp {input.genome} --gffFnp {input.gff} \
		--knownAminoAcidChangesFnp {input.targets_of_interest} \
		--dout {params.output_folder} --overWriteDir  --metaFnp meta.tsv \
		--numThreads {threads} --variantOccurrenceCutOff 1 --getPairwiseComps
		'''

rule rename_genes:
	'''
	converts gene IDs to gene names for more intuitive outputs
	'''
	input:
		called_aa_variants=config['output_folder']+'/variantCallOnSeqAndProtein/reports/summary/AAChangesInfo.tsv.gz',
		geneid_to_genename=config['geneid_to_genename']
	output:
		renamed_aa_variants=config['output_folder']+'/renamed_AAChangesInfo.tsv'
	script:
		'scripts/rename_genes.py'

rule variants_to_prevalences:
	'''
	converts variants to prevalences
	'''
	input:
		renamed_aa_variants=config['output_folder']+'/renamed_AAChangesInfo.tsv',
		metadata_sheet=config['metadata_sheet']
	params:
		coverage_threshold=config['coverage_threshold'],
		alternate_threshold=config['alternate_threshold'],
		summarize_by=config['summarize_by'],
		sample_column=config['sample_column'],
		latitude_name=config['latitude_name'],
		longitude_name=config['longitude_name'],
		output_dir=config['output_folder']
	output:
		prevalence_table=config['output_folder']+'/amino_acid_prevalences.tsv'
	script:
		'scripts/calculate_prevalences.py'