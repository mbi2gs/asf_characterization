## get_phylum_nog_dicts.py
# write json files containing dictionaries with tax_id:nog_list mappings

import sys
import csv
import json
import pandas as pd
csv.field_size_limit(sys.maxsize) # Required to load member files

# Converts a dict of NOG:[tax_ids] to tax_id:[NOGs]
def flip_member_dict(member_dict):
	NOGs_in_species = {}
	for key in member_dict.keys():
		tax_ids = member_dict[key]
		for tax_id in tax_ids:
			if tax_id in NOGs_in_species.keys():
				NOGs_in_species[tax_id].append(key) # add NOG to list of species' NOGS
			else:
				NOGs_in_species[tax_id] = [key]
	return NOGs_in_species

# Create a dict of NOG:[tax_ids] from the .members file available from eggNOG
def get_nog_members(nog_member_file):
	NOG_members = {} # dict of NOG_IDS:[species_tax_ids]
	with open(nog_member_file) as tsvfile:
		members = csv.reader(tsvfile, delimiter='\t')
		for line in members:
			nog = line[1]
			proteins = line[5].split(',')
			tax_ids = []
			for protein in proteins:
				tax_ids.append(protein.split('.')[0])
			NOG_members[nog] = tax_ids
	return NOG_members

# Get all NOG annotations in bactNOG
nog_annotations = pd.read_csv('../data/bactNOG/bactNOG.annotations.tsv',\
                              sep='\t',header=None,\
                              names=['bactNOG','group_name','protein_count',\
                                    'species_count','functional_category',\
                                    'functional_description'])
nog_annotations = nog_annotations.drop(nog_annotations.columns[[0]], axis=1)

# Read file describing number of metagenome samples that NOG is present in
nog_sample_count = pd.read_csv('../results/metagenome_samples_per_NOG.tsv',\
                              sep='\t',header=None,skiprows=1,\
                              names=['group_name','metagenome_sample_count'])
nog_annotations_with_count = pd.merge(nog_annotations,nog_sample_count,on='group_name')

mg_nogs = nog_annotations_with_count[nog_annotations_with_count['functional_category'].map(len) == 1]

# Get all NOGs and the species they belong to for each phylum (bacteroidetes,
# firmicutes, deferibacteres)
bcto_NOG_members = get_nog_members('../data/bctoNOG/bctoNOG.members.tsv')
NOGs_in_bcto = flip_member_dict(bcto_NOG_members)
firm_NOG_members = get_nog_members('../data/firmNOG/firmNOG.members.tsv')
NOGs_in_firm = flip_member_dict(firm_NOG_members)
# Get bact_nog_members to get list of all NOGs
bact_NOG_members = get_nog_members('../data/bactNOG/bactNOG.members.tsv')


# initialize dicts for tax_id:[NOGs] for each phylum
bcto_nogs_dict = {}
firm_nogs_dict = {}
def_nogs_dict = {}

# if a bactNOG is in any metagenome, check for it in every species
# if the NOG is present in a species within a phylum, add it to that tax_id entry
for NOG in bact_NOG_members.keys():
    if any(mg_nogs['group_name'] == NOG):
        for tax_id in NOGs_in_bcto.keys():
    		if tax_id in bact_NOG_members[NOG]:
    			if tax_id in bcto_nogs_dict.keys():
    				bcto_nogs_dict[tax_id].append(NOG) # dict of tax_id:bactNOGs
    			else:
    				bcto_nogs_dict[tax_id] = [NOG]
    	for tax_id in NOGs_in_firm.keys():
    		if tax_id in bact_NOG_members[NOG]:
    			if tax_id in firm_nogs_dict.keys():
    				firm_nogs_dict[tax_id].append(NOG)
    			else:
    				firm_nogs_dict[tax_id] = [NOG]


# Save phylum dictionaries as json files
# json.dump(bact_nogs_dict, open('../data/phylum_nogs/bact_nogs_dict.json','w'))
json.dump(bcto_nogs_dict, open('../data/phylum_nogs/bcto_nogs_dict.json','w'))
json.dump(firm_nogs_dict, open('../data/phylum_nogs/firm_nogs_dict.json','w'))
