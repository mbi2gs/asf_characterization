# asf_coverage_per_sample.py
# Calculates the coverage of metagenomic NOGs by the ASF per sample, per-category.
from __future__ import division

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pylab import rcParams
rcParams['figure.figsize'] = 10, 5 # default figure size

# Load metagenome data and NOG metadata
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

mg_nogs = nog_annotations_with_count[nog_annotations_with_count['functional_category'].map(len) == 1]


# Load ASF NOG coverage
asf_cover = pd.read_csv('../results/metagenome_coverage/asf_mg_coverage.tsv',sep='\t')

os.chdir('../results/metagenome_NOG_lists')
mg_nogs_by_sample = {}
# For each sample, get list of all nogs
for file in os.listdir(os.getcwd()):
    fname = file.split('/')[-1].split('.txt')[0]
    mg_nogs_by_sample[fname] = pd.read_csv(file,header=None,names=['group_name'])
    mg_nogs_by_sample[fname] = pd.merge(mg_nogs_by_sample[fname],mg_nogs,on="group_name")

# walk back to bin
os.chdir('../../bin')

# For each sample, calculate coverage in each category
ordered_categories = ['J','K','L','D','V','T','M','N','U','O','C'\
                                ,'G','E','F','H','I','P','Q','S']
sample_cover = []
for sample in mg_nogs_by_sample:
    # remove ASF coverage columns that represent nogs not present in sample
    coverage = {}
    for i,category in enumerate(ordered_categories):
        nogs_to_keep = mg_nogs_by_sample[sample][mg_nogs_by_sample[sample]['functional_category'] == category]
        mg_nogs_in_category = nogs_to_keep['group_name']
        asf_category_cover = asf_cover[list(set(mg_nogs_in_category) & set(asf_cover.columns.values.tolist()))]
        asf_category_cover = asf_category_cover.sum(axis=0)
        total_nogs_possible = len(mg_nogs_in_category)
        nogs_covered_by_asf = len(asf_category_cover)
        coverage[category] = nogs_covered_by_asf/total_nogs_possible
    sample_cover.append(coverage)
print sample_cover

# flip sample cover so that each array in the list is a category (with 15 data points each)
# number of elements in list should equal the number of categories
cover_by_category = {}
for category in ordered_categories:
    cover_by_category[category] = np.empty(15)
    for i,sample in enumerate(sample_cover):
        cover_by_category[category][i] = sample[category]

# Convert from dict to list and list of arrays and maintain pairing
category_list = []
value_list = []
for category in ordered_categories:
    category_list.append(category)
    value_list.append(cover_by_category[category])

# Create boxplot
numboxes = 19
plt.figure()
boxprops = dict(color='black',alpha=0.2)
whiskerprops = dict(linestyle='-',color='black',alpha = 0.2)
capprops = dict(color='black',alpha = 0.2)
medianprops = dict(color='blue',alpha= 0.5)
boxplot = plt.boxplot(value_list,labels=category_list,sym='',whis=[5,95],\
                        boxprops=boxprops,whiskerprops=whiskerprops,\
                        capprops=capprops,medianprops=medianprops,\
                        patch_artist=True)
# Set box color to black
for patch in boxplot['boxes']:
    patch.set_facecolor('black')
# Add scatter points
for i in range(numboxes):
    y = value_list[i]
    x = np.random.normal(1+i, 0.04, size=len(y))
    plt.plot(x,y,'r.',alpha=0.4)
plt.title('Coverage of metagenomic samples by ASF')
plt.xlabel('Functional Category')
plt.ylabel('Fraction of NOGs represented by ASF')
plt.savefig('../results/mg_cover_by_asf.svg',type='svg')
plt.show()
