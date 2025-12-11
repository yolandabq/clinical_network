#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 16:41:16 2025

@author: yolanda
"""

import pandas as pd
from itertools import combinations
from scipy.stats import zscore
import math
import os
import glob

#### functions

def jaccard_similarity(set1, set2):
    # Efficient Jaccard similarity calculation
    intersection = len(set1 & set2)
    union = len(set1 | set2)
    return intersection / union if union != 0 else 0

###########

prefix_folder = '/home/yolanda/'
#prefix_folder = '/mnt/'

panel_app_folder = str(prefix_folder + 'tblab/yolanda/GLOWgenes/berta/networks/PanelApp/')
panel_app_list = glob.glob(str(panel_app_folder + "*GAR.csv"))

df_panels = pd.DataFrame(columns=['gene', 'entity_type', 'panel_name'])

for panel_id in panel_app_list:
    df_panel = pd.read_table(panel_id)
    df_panel["panel_name"] = os.path.basename(panel_id).replace('_GAR.csv', '')
    df_panels = pd.concat([df_panels, df_panel.loc[:,["gene", "entity_type", "panel_name"]]])


print("Number of panels: " + str(len(set(df_panels["panel_name"]))))
print("Number of unique genes: " + str(len(set(df_panels["gene"]))))

gene_panel_dict = {}

for gen in set(df_panels["gene"]):
    gene_panel_dict[gen] = list(df_panels[df_panels["gene"] == gen]["panel_name"])

##### calculate jaccard between all pairs

similarities = {}

for gene1, gene2 in combinations(gene_panel_dict.keys(), 2):
    sim = jaccard_similarity(set(gene_panel_dict[gene1]), set(gene_panel_dict[gene2]))
    similarities[(gene1, gene2)] = sim


# Convert to a DataFrame
df_jaccard_sim = pd.DataFrame(
    [(g1, g2, sim) for (g1, g2), sim in similarities.items()],
    columns=['Gene1', 'Gene2', 'JaccardSimilarity']
)

####### calculate zscore

df_jaccard_sim = df_jaccard_sim[df_jaccard_sim['JaccardSimilarity']>0] # remove no jaccard 

df_jaccard_sim['Jaccard_Z'] = zscore(df_jaccard_sim['JaccardSimilarity'])

df_jaccard_sim.to_csv(str(prefix_folder + 'tblab/yolanda/GLOWgenes/berta/networks/PanelApp/PanelApp_jaccard.csv'), index=False)


############### filter significant edges and create cytoscape network

significant_edges = df_jaccard_sim[df_jaccard_sim['Jaccard_Z'] > 1.96]
print(f"Found {len(significant_edges)} significant gene-gene interactions")

cytoscape_edges = significant_edges[['Gene1', 'Gene2', 'Jaccard_Z']].copy()
cytoscape_edges.rename(columns={
    'Gene1': 'Source',
    'Gene2': 'Target',
    'Jaccard_Z': 'Weight'
}, inplace=True)


cytoscape_edges.to_csv(str(prefix_folder + 'tblab/yolanda/GLOWgenes/berta/networks/PanelApp/PanelApp_network_1_96.csv'), index=False)



