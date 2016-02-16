# principle_coordinate_analysis.py
# Script for generating PCoA plots

import pandas as pd
import matplotlib.pyplot as plt
import json
from sklearn.metrics import jaccard_similarity_score
from sklearn.metrics.pairwise import pairwise_distances
from skbio.stats.ordination import pcoa

# Compute Jaccard distances
# bcto_matrix = bcto_cover.as_matrix()
# bcto_distances = pairwise_distances(bcto_matrix,metric='jaccard')
# pcoa_of_distance = pcoa(bcto_distances)
# pcoa_of_distance.plot()
# plt.show()
# print pcoa_of_distance

# firm_matrix = firm_cover.as_matrix()
# firm_distances = pairwise_distances(firm_matrix,metric='jaccard')
# pcoa_of_distance = pcoa(firm_distances)
# pcoa_of_distance.plot()
# plt.show()
# print pcoa_of_distance

asf_bcto_firm_cover = pd.read_csv('../results/metagenome_coverage/asf_bcto_firm_mg_coverage.tsv',sep='\t')


asf_bcto_firm_matrix = asf_bcto_firm_cover.as_matrix()
asf_bcto_firm_distances = pairwise_distances(asf_bcto_firm_matrix,metric='jaccard')
pcoa_of_distance = pcoa(asf_bcto_firm_distances)
pcoa_of_distance.plot()
plt.show()
