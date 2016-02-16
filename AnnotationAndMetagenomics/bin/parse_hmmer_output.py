# parse_asf_hmmer_output.py
# script to parse table output from hmmer to get a list of all best-hit NOGs
# that occur at least once for each species.

from __future__ import division # forces int division to be flo
import os
import csv
import sys
from collections import Counter
import copy
import random
import pandas as pd

def get_nog_counts(directory,merge_tables=True):
	'''
	Returns dict of NOG counts if merge_tables=True, else returns dict of dict of
	NOG counts. The first situation is for merging multiple tables into a single
	result (ala metagenome files), whereas the second is for reporting the results
	from multiple single files (ala each file is a species)
	'''
	master_dict = {}
	master_dict = Counter(master_dict)
	os.chdir(directory)
	for table_file in os.listdir(directory):
		print table_file
		table = csv.reader(open(table_file,'r'), delimiter='\t')
		output_dict = {}
		max_score = 0
		min_e_val = 10
		best_hit = ''
		current_gene = ''
		line_count = 0
		for line in table:
			split_line = line[0].split()
			if split_line[0].find('#') < 0:
				if len(split_line[0].split('.')) < 2:
					print split_line[0].split('.')
				hit = split_line[0].split('.')[1] #Grab only the ENOG
				gene = split_line[2]
				e_val = float(split_line[4])
				score = float(split_line[5])
				if gene != current_gene:
					if len(current_gene) > 1 and min_e_val < 1E-10: #if this isn't the first gene, then a
											  #new gene indicates there are no more
											  #entries for the old gene.
						output_dict[best_hit] = 1
					current_gene = gene
					max_score = 0
					best_hit = hit
					new_gene = True
				if e_val < min_e_val and e_val < 1E-10:
					min_e_val = e_val
					best_hit = hit

				line_count += 1
			else:
				line_count += 1
		if merge_tables:
			master_dict = master_dict + Counter(output_dict)
		else:
			master_dict[table_file] = output_dict
	return master_dict


if __name__ == '__main__':
	os.chdir('../data/asf_hmmer_output')
	asf_nogs = get_nog_counts(os.getcwd(),merge_tables=False)
	os.chdir('../../bin')
	### write asf_nogs data.
	# output is a dict of dicts. Each dict is a species, therefore write
	# a tsv for each dict
	for species in asf_nogs.keys():
		with open('../results/asf_NOG_lists/all_nogs/' + species + '_nogs.txt','w') as f:
			for key in asf_nogs[species]:
				print >>f, key

	# Get NOGs from metagenomes
	os.chdir('../data/metagenome_hmmer_output/')
	sub_dirs = [x[0] for x in os.walk(os.getcwd())]
	sub_dirs = sub_dirs[1:] #remove the parent directory
	num_samples_NOG_is_in = {}
	for directory in sub_dirs:
		print directory
		print os.getcwd()
		sample_name = directory.split('/')[-1]
		print directory, sample_name
		NOGs_in_metagenome = get_nog_counts(directory).keys() # extract list of NOGS from keys of (dict of NOG:# of occurences).
		os.chdir('..')

		for NOG in NOGs_in_metagenome:
			if NOG in num_samples_NOG_is_in.keys():
				num_samples_NOG_is_in[NOG] += 1
			else:
				num_samples_NOG_is_in[NOG] = 1
		print os.getcwd()
		with open('../../results/metagenome_NOG_lists/' + sample_name + '_nogs.txt','w') as f:
			for entry in NOGs_in_metagenome:
				print >>f, entry

	metagenome_sample_count = pd.DataFrame.from_dict(num_samples_NOG_is_in,\
													 orient='index')
	metagenome_sample_count.to_csv('../../results/metagenome_samples_per_NOG.tsv',sep='\t')
else:
	print 'loading parse_hmmer_output.py'
