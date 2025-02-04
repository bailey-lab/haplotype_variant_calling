import os
import gzip
import subprocess

seekdeep_folder=snakemake.input.seekdeep_folder
amplicons=os.listdir(seekdeep_folder)
amplicons.remove('locationByIndex')
output_path=snakemake.output.seekdeep_tables[:-3] #trim off the .gz at the end
output_file=open(output_path, 'w')

for amplicon_number, amplicon in enumerate(sorted(amplicons)):
	cluster_file=f'{seekdeep_folder}/{amplicon}/analysis/selectedClustersInfo.tab.txt.gz'
	for line_number, line in enumerate(gzip.open(cluster_file, mode='rt')):
		if line_number==0:
			print('number of columns in', amplicon, 'is', len(line.strip().split('\t')))
		if line_number==0 and amplicon_number==0:
			output_file.write(line)
		elif line_number>0:
			output_file.write(line)
output_file.close()
subprocess.call(['gzip', '-k', output_path])