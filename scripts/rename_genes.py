'''
renames geneID transcripts to human-readable gene names
'''

aa_variants=snakemake.input.called_aa_variants
conversion_table=snakemake.input.geneid_to_genename
output_file=open(snakemake.output.renamed_aa_variants, 'w')
import gzip

conversion_dict=dict([line.strip().split('\t') for line in open(conversion_table)][1:])

for line_number, line in enumerate(gzip.open(aa_variants, mode='rt')):
	line=line.strip().split('\t')
	if line_number==0:
		gene_col=line.index('Gene')
		mutation_col=line.index('Mutation_Name')
	else:
		gene_id=line[gene_col]
		if gene_id in conversion_dict:
			new_name=conversion_dict[gene_id]
			line[gene_col]=line[gene_col].replace(gene_id, new_name)
			line[mutation_col]=line[mutation_col].replace(gene_id, new_name)
	output_file.write('\t'.join(line)+'\n')