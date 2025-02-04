'''
takes a table of counts and converts it to prevalences using thresholds. This
code is modified from the baileylab calculate_prevalences github repo
'''
aa_variants=snakemake.input.renamed_aa_variants
metadata_sheet=snakemake.input.metadata_sheet
coverage_threshold=snakemake.params.coverage_threshold
alternate_threshold=snakemake.params.alternate_threshold
summarize_by=snakemake.params.summarize_by
sample_column=snakemake.params.sample_column
prevalence_table=snakemake.output.prevalence_table
latitude_name=snakemake.params.latitude_name
longitude_name=snakemake.params.longitude_name
output_dir=snakemake.params.output_dir

def make_metadata_dict(metadata_sheet, sample_column, summarize_by, latitude_name, longitude_name):
	'''
	parses through the metadata sheet and assigns samples to categories using
	the summarize_by column. Also calculates average latitude and longitude for
	each category
	'''
	metadata_dict, coord_dict={},{}
	for line_number, line in enumerate(open(metadata_sheet)):
		if len(line.strip().split('\t'))>1:
			line=line.strip().split('\t')
		elif '"' not in line and "'" not in line:
			line=line.strip().split(',') #hopefully there are no commas in the metadata cells
		else:
			print('this program will not work if there are commas or quotes in the metadata cells. Check these carefully')
			exit()
		if line_number==0:
			sample_column=line.index(sample_column)
			category_column=line.index(summarize_by)
			latitude_column=line.index(latitude_name)
			longitude_column=line.index(longitude_name)
		else:
			sample, category=line[sample_column], line[category_column]
			latitude, longitude=line[latitude_column], line[longitude_column]
			metadata_dict[sample]=category
			coord_dict.setdefault(category, [[],[]])
			coord_dict[category][0].append(float(latitude))
			coord_dict[category][1].append(float(longitude))
	for category in coord_dict:
		coord_dict[category][0]=round(sum(coord_dict[category][0])/len(coord_dict[category][0]), 5)
		coord_dict[category][1]=round(sum(coord_dict[category][1])/len(coord_dict[category][1]), 5)
	return metadata_dict, coord_dict

def summarize_prevalences(metadata_dict, coverage_threshold, alternate_threshold, aa_variants):
	'''
	summarizes prevalences by checking if alternate counts and coverage counts
	are above thresholds for each sample, and then checking the category that
	each samples belongs to, and adding one to the summary dict if the sample
	surpasses the thresholds in a given category. Output is the number of
	samples that had each mutation and number of samples that had enough
	coverage in each sample
	'''
	#TODO: add a new function for checking coverage levels using final counts
	summary_dict, mutation_set, count_dict={}, set([]), {}
	header_dict={}
	for line_number, line in enumerate(open(aa_variants)):
		line=line.strip().split('\t')
		if line_number==0:
			sample_col=line.index('sample')
			gene_col=line.index('Gene')
			AA_pos_col=line.index('reference_AA_pos')
			ref_cnt_col=line.index('reference_AA_cnt')
			alt_cnt_col=line.index('alternate_AA_cnt')
			cov_cnt_col=line.index('coverage_AA_cnt')
			mutation_col=line.index('Mutation_Name')
			AA_change_col=line.index('AA_Change')
			targeted_col=line.index('Targeted')
			geneid_col=line.index('Gene_ID')
			exonic_func_col=line.index('ExonicFunc')
		else:
			sample=line[sample_col]
			gene=line[gene_col]
			AA_pos=int(line[AA_pos_col])
			gene_pos=gene+'_'+str(AA_pos)
			ref_cnt=int(line[ref_cnt_col])
			exonic_func=line[exonic_func_col]
			alt_cnt=int(line[alt_cnt_col])
			cov_cnt=int(line[cov_cnt_col])
			targeted=line[targeted_col]
			geneid=line[geneid_col]
			mutation=line[AA_change_col]
			category=metadata_dict[sample]
			gene_mut=line[mutation_col]
			mutation_set.add((gene, AA_pos, gene_mut))
			header_dict[gene_mut]=[geneid, gene, gene_mut, exonic_func, mutation, targeted]
			summary_dict.setdefault(category, {})
			summary_dict[category].setdefault(gene_mut, [0,0])
			count_dict.setdefault(sample, {})
			count_dict[sample].setdefault(gene_mut, {'reference':0, 'alternate':0, 'coverage':0, 'other':0})
			count_dict[sample][gene_mut]['reference']+=(ref_cnt)
			count_dict[sample][gene_mut]['alternate']+=(alt_cnt)
			count_dict[sample][gene_mut]['coverage']+=(cov_cnt)
			count_dict[sample][gene_mut]['other']+=cov_cnt-(ref_cnt+alt_cnt)
			if cov_cnt>coverage_threshold:
				summary_dict[category][gene_mut][1]+=1
				if alt_cnt>alternate_threshold:
					summary_dict[category][gene_mut][0]+=1
	sorted_mutations=[item[-1] for item in sorted(list(mutation_set))]
	return summary_dict, count_dict, sorted_mutations, header_dict

def write_prevalences(summary_dict, prevalence_table, sorted_mutations):
	'''
	for each category, iterates through mutations and outputs the number of
	samples that had the mutation of interest and the number of samples that had
	enough coverage to make a call.
	'''
	categories=sorted(list(summary_dict.keys()))
	output_file=open(prevalence_table, 'w')
	output_file.write(summarize_by+'\t'+'\t'.join(sorted_mutations)+'\n')
	for category in categories:
		output_line=[category]
		for mutation in sorted_mutations:
			numerator, denominator=summary_dict[category][mutation]
			if denominator>0:
				string_version=f'{round(numerator/denominator, 3)} ({numerator}/{denominator})'
			else:
				string_version=f'0.000 ({numerator}/{denominator})'
			output_line.append(string_version)
		output_file.write('\t'.join(output_line)+'\n')

def make_counts_tables(count_dict, sorted_mutations, header_dict, output_dir):
	'''
	converts Nick's tables into Ozkan format for use with downstream parsers.
	'''
	file_types=['reference', 'alternate', 'coverage', 'other']
	output_files=[open(f'{output_dir}/{file_type}_AA_table.csv', 'w') for file_type in file_types]
	titles=['Gene ID', 'Gene', 'Mutation Name', 'ExonicFunc', 'AA Change', 'Targeted']
	for pos in range(6):
		output_line=[titles[pos]]
		for mutation in sorted_mutations:
			output_line.append(header_dict[mutation][pos])
		for output_file in output_files:
			output_file.write(','.join(output_line)+'\n')
	for sample in sorted(list(count_dict.keys())):
		for file_number in range(len(output_files)):
			output_line=[sample]
			for mutation in sorted_mutations:
				output_line.append(str(count_dict[sample][mutation][file_types[file_number]]))
			output_files[file_number].write(','.join(output_line)+'\n')

metadata_dict, coord_dict=make_metadata_dict(metadata_sheet, sample_column, summarize_by, latitude_name, longitude_name)
summary_dict, count_dict, sorted_mutations, header_dict=summarize_prevalences(metadata_dict, coverage_threshold, alternate_threshold, aa_variants)
write_prevalences(summary_dict, prevalence_table, sorted_mutations)
make_counts_tables(count_dict, sorted_mutations, header_dict, output_dir)