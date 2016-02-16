## asf_NOG_analysis.py
# analysis of NOG lists from each ASF species
from __future__ import division # forces int division to be flo
import os
import copy
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sortedcontainers import SortedSet
from pylab import rcParams
rcParams['figure.figsize'] = 10, 5 # default figure size

# Get bacterial NOG annotations and metagenome sample counts for each NOG
nog_annotations = pd.read_csv('../data/bactNOG/bactNOG.annotations.tsv',\
                              sep='\t',header=None,\
                              names=['bactNOG','group_name','protein_count',\
                                    'species_count','functional_category',\
                                    'functional_description'])
nog_annotations = nog_annotations.drop(nog_annotations.columns[[0]], axis=1)
nog_sample_count = pd.read_csv('../results/metagenome_samples_per_NOG.tsv',\
                              sep='\t',header=None,skiprows=1,\
                              names=['group_name','metagenome_sample_count'])
nog_annotations_with_count = pd.merge(nog_annotations,nog_sample_count,on='group_name')
# Parse each ASF NOG file and create a new dataframe for each species.
asf_nogs = {}
nog_counts = {}
for nog_file in os.listdir('../results/asf_NOG_lists/all_nogs/'):
    asf_nogs[nog_file[:4]] = pd.read_csv('../results/asf_NOG_lists/all_nogs/' + nog_file,header=None,names=['group_name'])
    asf_nogs[nog_file[:4]] = pd.merge(nog_annotations,asf_nogs[nog_file[:4]],on='group_name')
    for nog in asf_nogs[nog_file[:4]]['group_name']:
        if nog not in nog_counts:
            nog_counts[nog] = 1
        else:
            nog_counts[nog] += 1

# Get Unique NOGs
unique_nogs = {} # map nog to species it occurs in
for nog in nog_counts.keys():
    if nog_counts[nog] == 1: # if nog is unique, find species it's from
        for species in asf_nogs.keys():
            if nog in asf_nogs[species]['group_name'].values:
                unique_nogs[nog] = species

unique_nogs = pd.DataFrame({'group_name':unique_nogs.keys(),'asf_species':unique_nogs.values()})
unique_nogs_with_annotations = pd.merge(nog_annotations_with_count,unique_nogs,on='group_name')
unique_nogs_with_annotations.to_csv('../results/asf_unique_nogs.tsv', sep='\t')

unique_nogs_in_core = copy.deepcopy(unique_nogs_with_annotations)
unique_nogs_in_core = unique_nogs_in_core[unique_nogs_in_core['metagenome_sample_count'] == 15]
unique_nogs_in_core.to_csv('../results/asf_nogs_in_core_metagenome.tsv', sep='\t')



# create file for each ASF species containing NOGs that are present in metagenomes
asf_nog_dicts = {}
for species in asf_nogs:
    nogs = asf_nogs[species]
    asf_nog_dicts[species] = nogs['group_name'].tolist()
json.dump(asf_nog_dicts, open('../results/asf_NOG_lists/asf_nog_dicts.json','w'))
