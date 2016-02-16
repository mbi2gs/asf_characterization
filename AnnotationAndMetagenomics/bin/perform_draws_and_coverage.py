# jaccard.py
# Computes the coverage of all bactNOG species over the metagenome, then Computes
# the jaccard distance between every species and performs PCoA.
from __future__ import division
import pandas as pd
import matplotlib.pyplot as plt
import json
import numpy as np
from pylab import rcParams
rcParams['figure.figsize'] = 10, 5 # default figure size


# Get metagenome NOGs
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

# Load NOG coverage
asf_cover = pd.read_csv('../results/metagenome_coverage/asf_mg_coverage.tsv',sep='\t')
bcto_cover = pd.read_csv('../results/metagenome_coverage/bcto_mg_coverage.tsv',sep='\t')
firm_cover = pd.read_csv('../results/metagenome_coverage/firm_mg_coverage.tsv',sep='\t')

# Compute ASF np array for % coverage
# initialize np array before loop
asf_coverage = np.empty(15)
for i in range(1,16):
    nogs_w_sample_count = mg_nogs[mg_nogs['metagenome_sample_count'] == i]
    nogs_to_keep = nogs_w_sample_count['group_name'].tolist() #list of nogs that occur in i # of metagenomes
    # Get asf coverage of nogs
    asf_sub_cover = asf_cover[list(set(nogs_to_keep) & set(asf_cover.columns.values.tolist()))]
    asf_sub_cover = asf_sub_cover.sum(axis=0)
    # Get fraction of non-zero entries, then store in np array
    total_nogs_possible = len(nogs_to_keep)
    nogs_covered_by_asf = len(asf_sub_cover)

    print "total possible: " + str(total_nogs_possible),"covered: " + str(nogs_covered_by_asf)
    asf_coverage[i-1] = nogs_covered_by_asf/total_nogs_possible


num_draws = 10000
draw_coverage = np.empty([num_draws,15])
large_draw_draw_coverage = np.empty([num_draws,15])
larger_draw_draw_coverage = np.empty([num_draws,15])

for draw in range(0,num_draws):
    bcto_species = bcto_cover.sample(n=2)
    firm_species = firm_cover.sample(n=6)
    large_draw_bcto_species = bcto_cover.sample(n=4)
    large_draw_firm_species = firm_cover.sample(n=12)
    larger_draw_bcto_species = bcto_cover.sample(n=8)
    larger_draw_firm_species = firm_cover.sample(n=24)

    # Rmoved "Unnamed: 0" from parsed cover file
    bcto_species.drop('Unnamed: 0', axis=1, inplace=True)
    firm_species.drop('Unnamed: 0', axis=1, inplace=True)
    large_draw_bcto_species.drop('Unnamed: 0', axis=1, inplace=True)
    large_draw_firm_species.drop('Unnamed: 0', axis=1, inplace=True)
    larger_draw_bcto_species.drop('Unnamed: 0', axis=1, inplace=True)
    larger_draw_firm_species.drop('Unnamed: 0', axis=1, inplace=True)

    # drop NOG columns in species dataframes which have all 0's (i.e. nog is not present in any of the selected species)
    bcto_species = bcto_species.loc[:,(bcto_species != 0).any(axis=0)]
    firm_species = firm_species.loc[:,(firm_species != 0).any(axis=0)]
    large_draw_bcto_species = large_draw_bcto_species.loc[:,(large_draw_bcto_species != 0).any(axis=0)]
    large_draw_firm_species = large_draw_firm_species.loc[:,(large_draw_firm_species != 0).any(axis=0)]
    larger_draw_bcto_species = larger_draw_bcto_species.loc[:,(larger_draw_bcto_species != 0).any(axis=0)]
    larger_draw_firm_species = larger_draw_firm_species.loc[:,(larger_draw_firm_species != 0).any(axis=0)]
    for i in range(1,16):
        nogs_w_sample_count = mg_nogs[mg_nogs['metagenome_sample_count'] == i]
        nogs_to_keep = nogs_w_sample_count['group_name'].tolist() #list of nogs that occur in i # of metagenomes
        total_nogs_possible = len(nogs_to_keep)

        # Get coverage: intersection of the list of NOGs from metagenomes and draws
        bcto_sub_cover = bcto_species[list(set(nogs_to_keep) & set(bcto_species.columns.values.tolist()))]
        firm_sub_cover = firm_species[list(set(nogs_to_keep) & set(firm_species.columns.values.tolist()))]
        draw_cover = list(set(firm_sub_cover.columns.values.tolist()) | set(bcto_sub_cover.columns.values.tolist()))

        nogs_covered_by_draw = len(draw_cover)
        draw_coverage[draw,i-1] = nogs_covered_by_draw/total_nogs_possible

        # Get coverage for large draw set
        large_draw_bcto_sub_cover = large_draw_bcto_species[list(set(nogs_to_keep) & set(large_draw_bcto_species.columns.values.tolist()))]
        large_draw_firm_sub_cover = large_draw_firm_species[list(set(nogs_to_keep) & set(large_draw_firm_species.columns.values.tolist()))]
        large_draw_draw_cover = list(set(large_draw_firm_sub_cover.columns.values.tolist()) | set(large_draw_bcto_sub_cover.columns.values.tolist()))
        nogs_covered_by_draw = len(large_draw_draw_cover)
        large_draw_draw_coverage[draw,i-1] = nogs_covered_by_draw/total_nogs_possible

        larger_draw_bcto_sub_cover = larger_draw_bcto_species[list(set(nogs_to_keep) & set(larger_draw_bcto_species.columns.values.tolist()))]
        larger_draw_firm_sub_cover = larger_draw_firm_species[list(set(nogs_to_keep) & set(larger_draw_firm_species.columns.values.tolist()))]
        larger_draw_draw_cover = list(set(larger_draw_firm_sub_cover.columns.values.tolist()) | set(larger_draw_bcto_sub_cover.columns.values.tolist()))
        nogs_covered_by_draw = len(larger_draw_draw_cover)
        larger_draw_draw_coverage[draw,i-1] = nogs_covered_by_draw/total_nogs_possible

# calculate median coverage at each draw level
median_coverage = np.median(draw_coverage, axis=0)
large_draw_median_coverage = np.median(large_draw_draw_coverage, axis=0)
larger_draw_median_coverage = np.median(larger_draw_draw_coverage, axis=0)
# calculate 5th and 95th percentiles of coverage
fifth_percentile = np.percentile(draw_coverage,5,axis=0)
ninety_fifth_percentile = np.percentile(draw_coverage,95,axis=0)
large_draw_fifth_percentile = np.percentile(large_draw_draw_coverage,5,axis=0)
large_draw_ninety_fifth_percentile = np.percentile(large_draw_draw_coverage,95,axis=0)
larger_draw_fifth_percentile = np.percentile(larger_draw_draw_coverage,5,axis=0)
larger_draw_ninety_fifth_percentile = np.percentile(larger_draw_draw_coverage,95,axis=0)

# create array for x values (presence in 1:15 metagenomic samples)
sample_frequency = np.arange(1,16,1)
# fig, (ax,ax2) = plt.subplots(1,2,sharey=True)
fig, ax = plt.subplots(1,1,sharey=True)
ax.plot(sample_frequency,asf_coverage,label='ASF (8 species)',linewidth=1.0)
ax.plot(sample_frequency,median_coverage,linewidth=0.5,color='black',alpha=0.5)
ax.plot(sample_frequency,large_draw_median_coverage,linewidth=0.5,color='black',alpha=0.5)
ax.plot(sample_frequency,larger_draw_median_coverage,linewidth=0.5,color='black',alpha=0.5)
#  Adjust median line width
# ax2.plot(sample_frequency,asf_coverage,sample_frequency,large_draw_median_coverage)

ax.fill_between(sample_frequency,fifth_percentile,ninety_fifth_percentile,label='8 species',facecolor='red',alpha=0.4)
ax.fill_between(sample_frequency,large_draw_fifth_percentile,large_draw_ninety_fifth_percentile,label='16 species',facecolor='red',alpha=0.2)
ax.fill_between(sample_frequency,larger_draw_fifth_percentile,larger_draw_ninety_fifth_percentile,label='32 species',facecolor='red',alpha=0.05)
ax.legend(loc='upper left')
# ax2.fill_between(sample_frequency,large_draw_fifth_percentile,large_draw_ninety_fifth_percentile,facecolor='green',alpha=0.33)

ax.set_title('Metagenomic coverage: ASF vs. random draws')
ax.set_ylabel('Fraction of NOGs represented')
ax.set_xlabel('Frequency (# metagenomes NOG is present in)')
ax.set_xlim([1,15])

# ax2.set_title('Metagenomic coverage: ASF vs. random draws (80 species)')
# ax2.set_ylabel('Percentage of NOGs Covered')
# ax2.set_xlabel('Frequency (# metagenomes NOG is present in)')
# ax2.set_xlim([1,15])

plt.savefig('../results/asf_v_draws.svg',type='svg')
plt.show()
