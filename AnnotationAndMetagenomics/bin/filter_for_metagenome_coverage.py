# filter_for_metagenome_coverage.py
# Takes HMMer output and filters the output NOGs so that only those contained
# in at least one metagenome are returned.

import pandas as pd
import matplotlib.pyplot as plt
import json
from sklearn.metrics import jaccard_similarity_score
from sklearn.metrics.pairwise import pairwise_distances
from skbio.stats.ordination import pcoa

def filter_species_nog_dict(nog_dict,filter_by):
    filter_by = set(filter_by.tolist())
    for species in nog_dict:
        nog_dict[species] = set(nog_dict[species]).intersection(filter_by)
        nogs = nog_dict[species]
    return nog_dict

def convert_nog_dict_to_coverage(nog_dict):
    #convert list to dict of entry:1
    for species in nog_dict:
        nog_dict[species] = [{nog:1} for nog in nog_dict[species]]

    mg_nog_cover = {}
    for species in nog_dict:
        species_nog_cover = pd.DataFrame.from_dict(nog_dict[species])
        mg_nog_cover[species] = species_nog_cover.sum(axis=0)
    mg_nog_cover_df = pd.DataFrame.from_dict(mg_nog_cover)
    mg_nog_cover_df = mg_nog_cover_df.transpose()
    mg_nog_cover_df = mg_nog_cover_df.fillna(0)
    return mg_nog_cover_df

# open metagenome NOG file
mg_nog_annotations = pd.read_csv('../data/bactNOG/bactNOG.annotations.tsv',\
                              sep='\t',header=None,\
                              names=['bactNOG','group_name','protein_count',\
                                    'species_count','functional_category',\
                                    'functional_description'])
mg_nog_annotations = mg_nog_annotations.drop(mg_nog_annotations.columns[[0]], axis=1)
mg_nog_sample_count = pd.read_csv('../results/metagenome_samples_per_NOG.tsv',\
                              sep='\t',header=None,skiprows=1,\
                              names=['group_name','metagenome_sample_count'])
mg_nog_annotations_with_count = pd.merge(mg_nog_annotations,mg_nog_sample_count,on='group_name')

# Get NOG information for all bacteroidetes and firmicutes
bcto_nogs_dict = json.load(open('../data/phylum_nogs/bcto_nogs_dict.json'))
firm_nogs_dict = json.load(open('../data/phylum_nogs/firm_nogs_dict.json'))
asf_nogs_dict = json.load(open('../results/asf_NOG_lists/asf_nog_dicts.json'))

# Filter dictionaries to only include NOGs that are in metagenomic samples
bcto_nogs_dict = filter_species_nog_dict(bcto_nogs_dict,mg_nog_annotations['group_name'])
firm_nogs_dict = filter_species_nog_dict(firm_nogs_dict,mg_nog_annotations['group_name'])
asf_nogs_dict = filter_species_nog_dict(asf_nogs_dict,mg_nog_annotations['group_name'])

# Generate binary coverage dataframes
asf_bcto_firm_nogs_dict = bcto_nogs_dict.copy()
asf_bcto_firm_nogs_dict.update(firm_nogs_dict.copy())
asf_bcto_firm_nogs_dict.update(asf_nogs_dict.copy())
asf_cover = convert_nog_dict_to_coverage(asf_nogs_dict)
bcto_cover = convert_nog_dict_to_coverage(bcto_nogs_dict)
firm_cover = convert_nog_dict_to_coverage(firm_nogs_dict)
asf_bcto_firm_cover = convert_nog_dict_to_coverage(asf_bcto_firm_nogs_dict)

# Save binary coverage dataframes
asf_cover.to_csv('../results/metagenome_coverage/asf_mg_coverage.tsv',sep='\t')
bcto_cover.to_csv('../results/metagenome_coverage/bcto_mg_coverage.tsv',sep='\t')
firm_cover.to_csv('../results/metagenome_coverage/firm_mg_coverage.tsv',sep='\t')
asf_bcto_firm_cover.to_csv('../results/metagenome_coverage/asf_bcto_firm_mg_coverage.tsv',sep='\t')
