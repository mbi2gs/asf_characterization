## asf_NOG_analysis.py
# analysis of NOG lists from each ASF species
from __future__ import division # forces int division to be flo
import os
import copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sortedcontainers import SortedSet
from pylab import rcParams
rcParams['figure.figsize'] = 12, 5 # default figure size

def generate_figure(fig,nog_df,data_filename,figure_filename,ylabel,title,drop=[]):
    unique_nogs = copy.deepcopy(nog_df)
    functional_categories = list(SortedSet(unique_nogs['functional_category']))
    functional_categories = [entry for entry in functional_categories if len(entry) == 1]
    asf_species = ['A356','A360','A361','A457','A492','A500','A502','A519']
    actual_names = ['ASF356','ASF360','ASF361','ASF457','ASF492','ASF500','ASF502','ASF519']
    convert_species = {'ASF356':'A356','ASF360':'A360','ASF361':'A361','ASF457':'A457','ASF492':'A492','ASF500':'A500','ASF502':'A502','ASF519':'A519'}
    df_by_category = pd.DataFrame(index=actual_names)
    # Fill matrix of asf_species rows, cateogry columns, with NOG count in entries
    for functional_category in functional_categories:
        df_by_category[functional_category] = 0
        for full_species in actual_names:
            species = convert_species[full_species]
            df_by_category[functional_category][full_species] = len(unique_nogs.loc[(unique_nogs['asf_species'] == species) & (unique_nogs['functional_category'] == functional_category)].index)
    if drop:
        for category in drop:
            df_by_category = df_by_category.drop(category,1) # remove category A since there is only 1 NOG
    df_by_category = df_by_category.reindex(index=actual_names)

    categories_by_class_without_missing = ['J','K','L','D','V','T','M','N','U','O','C'\
     								,'G','E','F','H','I','P','Q','S']
    tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),\
			(44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),\
			(148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),\
			(227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),\
			(188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
    for i in range(len(tableau20)):
		r, g, b = tableau20[i]
		tableau20[i] = (r/255, g/255, b/255)
    colors = {species:tableau20[i] for i,species in enumerate(actual_names)}
    df_by_category = df_by_category.transpose()
    df_by_category = df_by_category.reindex(index=categories_by_class_without_missing)
    df_by_category.to_csv('../results/' + data_filename + 'unnormed.tsv', sep='\t')
    df_by_category = df_by_category.div(df_by_category.sum(axis=1), axis=0)
    ax = df_by_category.plot(kind='bar', stacked=True, color=colors.values(), width=.85)
    df_by_category.to_csv('../results/' + data_filename + '.tsv', sep='\t')
    handles, labels = ax.get_legend_handles_labels()
    new_handles = []
    for i,handle in enumerate(handles):
        plt.setp(handle,fc=colors[labels[i]])
        new_handles.append(handle)
    ax.legend(handles=new_handles,labels=labels,loc='upper left')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1],bbox_to_anchor=(1.00, 0, 1., .101), loc=3, prop={'size':10})
    plt.xlabel('Functional Category')
    plt.ylabel(ylabel)
    plt.xticks(np.arange(len(categories_by_class_without_missing)), categories_by_class_without_missing, rotation='horizontal')
    plt.ylim(0,1)
    ax.set_title(title)
    ax.grid(False)

    plt.savefig('../results/' + figure_filename,type='svg', bbox_inches='tight')
    return ax

unique_nogs_with_annotations = pd.DataFrame.from_csv('../results/asf_unique_nogs.tsv', sep='\t')

fig = plt.figure()
ax1 = generate_figure(fig,unique_nogs_with_annotations,'fig6_data','figure_6.svg','Fraction of unique NOGs contributed','Unique NOGs in each ASF species',drop=[])

unique_nogs_in_core = pd.DataFrame.from_csv('../results/asf_nogs_in_core_metagenome.tsv', sep='\t')
ax2 = generate_figure(fig,unique_nogs_in_core,'fig7_data','figure_7.svg','Fraction of unique core NOGs contributed','Core NOGs contributed by each species')
plt.show()
