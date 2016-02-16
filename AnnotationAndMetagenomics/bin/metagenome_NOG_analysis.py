## metagenome_NOG_analysis.py
# analysis of the distribution of metagenome NOGs

import os
import time
import copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sortedcontainers import SortedSet

def generate_figure_1(fig,nog_annotations_with_count):
    df_by_known_unknown = []
    unknown_annotations = nog_annotations_with_count[nog_annotations_with_count['functional_category'] == 'S']
    known_annotations = nog_annotations_with_count[nog_annotations_with_count['functional_category'] != 'S']
    known_annotations.loc[:,'functional_category'] = 'Known'
    unknown_annotations.loc[:,'functional_category'] = 'Unknown'
    df_by_known_unknown.append(known_annotations)
    df_by_known_unknown.append(unknown_annotations)


    ax = fig.add_subplot(111)
    colors = ['b','r']
    n, bins, patches = ax.hist([df['metagenome_sample_count'] for df in df_by_known_unknown], \
        bins = np.arange(16)+0.5, stacked=True, color=colors, align='left')

    known_patch = mpatches.Patch(color=colors[0], label='Known function')
    unknown_patch = mpatches.Patch(color=colors[1], label='Unknown function')
    ax.legend(handles=[known_patch,unknown_patch], loc='upper right')
    plt.xlabel('Number of metagenomes NOG is present in')
    plt.ylabel('NOG count')
    plt.xticks(np.arange(15)+0.5, np.arange(15)+1, rotation='horizontal')
    ax.grid(False)
    plt.savefig('../results/figure_S1.eps',type='eps')
    return ax

def generate_figure_2(fig,nog_annotations_with_count,data_filename_1='fig1_data',data_filename_2='fig2_data'):
    unknown_annotations = nog_annotations_with_count[nog_annotations_with_count['functional_category'] == 'S']
    known_annotations = nog_annotations_with_count[nog_annotations_with_count['functional_category'] != 'S']
    known_annotations.loc[:,'functional_category'] = 'Known'
    unknown_annotations.loc[:,'functional_category'] = 'Unknown'

    norm_known_unknown = pd.DataFrame({\
        'Known':[len(known_annotations[known_annotations['metagenome_sample_count'] == i]) \
        for i in range(1,16)],\
        'Unknown':[len(unknown_annotations[unknown_annotations['metagenome_sample_count'] == i]) \
        for i in range(1,16)]})
    # normalize across each row to total 1
    norm_known_unknown.to_csv('../results/' + data_filename_1 + '.tsv', sep='\t')
    norm_known_unknown = norm_known_unknown.div(norm_known_unknown.sum(axis=1), axis=0)
    colors = ['b','r']
    ax = norm_known_unknown.plot(kind='bar', stacked=True, color=colors)
    norm_known_unknown.to_csv('../results/' + data_filename_2 + '.tsv', sep='\t')
    known_patch = mpatches.Patch(color=colors[0], label='Known function')
    unknown_patch = mpatches.Patch(color=colors[1], label='Unknown function')
    ax.legend(handles=[known_patch,unknown_patch], loc='upper right')
    plt.xlabel('Number of metagenomes NOG is present in')
    plt.ylabel('NOG abundance')
    plt.xticks(np.arange(15), np.arange(15)+1, rotation='horizontal')
    ax.grid(False)
    plt.savefig('../results/figure_S2.eps',type='eps')
    return ax

def generate_figure_3_4(fig,nog_annotations_with_count,data_filename_1='fig3_data',data_filename_2='fig4_data'):
    known_annotations = copy.deepcopy(nog_annotations_with_count)
    # replace S entries with 0 counts. Useful to keep S in to manage colors in legends
    known_annotations.loc[known_annotations['functional_category'] == 'S','metagenome_sample_count'] = 0
    functional_categories = list(SortedSet(known_annotations['functional_category']))
    functional_categories = [entry for entry in functional_categories if len(entry) == 1]
    df_by_category = []
    for functional_category in functional_categories:
        df_by_category.append(known_annotations[known_annotations['functional_category'] == functional_category])


    ax = fig.add_subplot(111)
    cm = plt.cm.get_cmap('Set3')
    colors = {category:cm(i*.04) for i,category in enumerate(functional_categories)}
    count_samples = pd.DataFrame({list(set(df['functional_category']))[0]:\
                    [len(df[df['metagenome_sample_count'] == i].index) for i in range(1,16)] for df in df_by_category})
    ax1 = count_samples.plot(kind='bar', stacked=True, color=colors.values())
    count_samples.to_csv('../results/' + data_filename_1 + '.tsv', sep='\t')
    handles, labels = ax1.get_legend_handles_labels()
    new_handles = []
    for i,handle in enumerate(handles):
        plt.setp(handle,fc=colors[labels[i]])
        new_handles.append(handle)
    ax1.legend(handles=new_handles,labels=labels,loc='upper right')
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles[::-1], labels[::-1],bbox_to_anchor=(1.00, 0, 1., .102), loc=3, prop={'size':10})
    plt.xlabel('Number of metagenomes NOG is present in')
    plt.ylabel('NOG count')
    plt.xticks(np.arange(15), np.arange(15)+1, rotation='horizontal')
    ax1.grid(False)
    plt.savefig('../results/figure_S3.eps',type='eps')

    count_samples = count_samples.div(count_samples.sum(axis=1), axis=0)
    ax2 = count_samples.plot(kind='bar', stacked=True, color=colors.values())
    count_samples.to_csv('../results/' + data_filename_2 + '.tsv', sep='\t')
    handles, labels = ax2.get_legend_handles_labels()
    new_handles = []
    for i,handle in enumerate(handles):
        plt.setp(handle,fc=colors[labels[i]])
        new_handles.append(handle)
    ax2.legend(handles=new_handles,labels=labels,loc='upper right')
    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles[::-1], labels[::-1],bbox_to_anchor=(1.00, 0, 1., .102), loc=3, prop={'size':10})

    plt.ylim(0,1)
    plt.xlabel('Number of metagenomes NOG is present in')
    plt.ylabel('NOG abundance')
    plt.xticks(np.arange(15), np.arange(15)+1, rotation='horizontal')
    ax2.grid(False)
    plt.savefig('../results/figure_S4.eps',type='eps')
    return ax1, ax2

def generate_figure_5(fig,nog_annotations_with_count,data_filename='fig5_data'):
    functional_categories = list(SortedSet(nog_annotations_with_count['functional_category']))
    functional_categories = [entry for entry in functional_categories if len(entry) == 1]
    core_nogs = nog_annotations_with_count[nog_annotations_with_count['metagenome_sample_count'] == 15]
    core_nogs_counts = {}
    all_nogs_counts = {}
    for functional_category in functional_categories:
        core_nogs_counts[functional_category] = [len(\
            core_nogs[core_nogs['functional_category'] == functional_category])]
        all_nogs_counts[functional_category] = [len(\
            nog_annotations_with_count[nog_annotations_with_count['functional_category'] == functional_category])]
    core_nogs_no_unknown_counts = copy.deepcopy(core_nogs_counts)
    all_nogs_no_unknown_counts = copy.deepcopy(all_nogs_counts)
    core_nogs_no_unknown_counts['S'] = 0
    all_nogs_no_unknown_counts['S'] = 0

    core_v_all = pd.DataFrame(all_nogs_counts)
    core_v_all = pd.concat([core_v_all,pd.DataFrame(core_nogs_counts)])
    core_v_all = pd.concat([core_v_all,pd.DataFrame(all_nogs_no_unknown_counts)])
    core_v_all = pd.concat([core_v_all,pd.DataFrame(core_nogs_no_unknown_counts)])
    core_v_all.index.names = ['NOG source']
    index = core_v_all.index
    names = index.names
    index = ['All','Core','All_no_unknown','Core_no_unknown']
    core_v_all.index = index
    core_v_all.fillna(0) # Unknown category will have NA's for _no_unknown_'s
    # print the number of total NOGs and the number of core NOGs.
    total_NOG_count = core_v_all.sum(axis=1)
    print '#####################################'
    print 'NOG counts:'
    print total_NOG_count
    norm_core_v_all = core_v_all.div(core_v_all.sum(axis=1),axis=0)
    cm = plt.cm.get_cmap('Set3')
    colors = {category:cm(i*.04) for i,category in enumerate(functional_categories)}
    core_v_all.to_csv('../results/' + data_filename + '.tsv', sep='\t')
    norm_core_v_all.to_csv('../results/' + data_filename + 'norm.tsv', sep='\t')
    ax = norm_core_v_all.plot(kind='bar',stacked=True, colors=colors.values())
    handles, labels = ax.get_legend_handles_labels()
    new_handles = []
    for i,handle in enumerate(handles):
        plt.setp(handle,fc=colors[labels[i]])
        new_handles.append(handle)
    ax.legend(handles=new_handles,labels=labels,loc='upper left')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1],bbox_to_anchor=(1.00, 0, 1., .102), loc=3, prop={'size':10})
    plt.xticks(np.arange(4), ['All','Core','All, unknown\nexcluded','Core, unknown\nexcluded'], rotation='horizontal')
    ax.grid(False)
    plt.ylim(0,1)
    plt.savefig('../results/figure_S5.eps',type='eps')
    return ax

# parse bactnog to get category for each NOG
nog_sample_count = pd.read_csv('../results/metagenome_samples_per_NOG.tsv',\
                              sep='\t',header=None,skiprows=1,\
                              names=['group_name','metagenome_sample_count'])
nog_annotations = pd.read_csv('../data/bactNOG/bactNOG.annotations.tsv',\
                              sep='\t',header=None,\
                              names=['bactNOG','group_name','protein_count',\
                                    'species_count','functional_category',\
                                    'functional_description'])
nog_annotations = nog_annotations.drop(nog_annotations.columns[[0]], axis=1)
nog_annotations_with_count = pd.merge(nog_annotations,nog_sample_count,on='group_name')


# Create df for each functional category
functional_categories = list(set(nog_annotations_with_count['functional_category']))
functional_categories = [entry for entry in functional_categories if len(entry) == 1]
df_by_category = []
for functional_category in functional_categories:
    df_by_category.append(nog_annotations_with_count[nog_annotations_with_count['functional_category'] == functional_category])

fig = plt.figure()
## Figure A: histogram of known vs. unknown functional annotations.
ax1 = generate_figure_1(fig,nog_annotations_with_count)
## Figure B: Figure A, normalized to 1 for all 15 x values
ax2 = generate_figure_2(fig,nog_annotations_with_count)
## Figure C & D: histogram of known NOGs only, abundances of known NOGs only
ax3,ax4 = generate_figure_3_4(fig,nog_annotations_with_count)
ax5 = generate_figure_5(fig,nog_annotations_with_count)
plt.show()
