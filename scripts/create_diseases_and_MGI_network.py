#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  4 10:49:11 2025

@author: yolanda
"""

import pandas as pd
from itertools import combinations
from scipy.stats import zscore
import math


#### functions

def jaccard_similarity(set1, set2):
    # Efficient Jaccard similarity calculation
    intersection = len(set1 & set2)
    union = len(set1 | set2)
    return intersection / union if union != 0 else 0

#################

url_diseases_annotation = "https://download.jensenlab.org/human_disease_textmining_filtered.tsv"
diseases_annotation = pd.read_table(url_diseases_annotation, header=None)

diseases_annotation.head()

#################

diseases_annotation = diseases_annotation.set_axis(['gene_identifier', 'gene_name', 'disease_identifier', 'disease_name', 'z-score', 'confidence_score', 'URL'], axis=1)
diseases_annotation.head()

## Filtramos? 
diseases_annotation_filt = diseases_annotation[diseases_annotation['z-score'] > 1.96]

################################
## numero de genes únicos
unique_genes = list(set(diseases_annotation_filt["gene_name"]))

print("Hay ", str(len(unique_genes)), " genes únicos.")

#######################################33

gene_disease_dict = {}

for gen in unique_genes:
    gene_disease_dict[gen] = list(diseases_annotation_filt[diseases_annotation_filt["gene_name"] == gen]["disease_identifier"])

print("En total, hay " + str(len(gene_disease_dict.keys())) + " genes únicos diferentes")
## En total, hay 20811 genes únicos diferentes

##### calculate jaccard between all pairs

similarities = {}

for gene1, gene2 in combinations(gene_disease_dict.keys(), 2):
    sim = jaccard_similarity(set(gene_disease_dict[gene1]), set(gene_disease_dict[gene2]))
    similarities[(gene1, gene2)] = sim


# Convert to a DataFrame
df_jaccard_sim = pd.DataFrame(
    [(g1, g2, sim) for (g1, g2), sim in similarities.items()],
    columns=['Gene1', 'Gene2', 'JaccardSimilarity']
)


####### calculate zscore

df_jaccard_sim = df_jaccard_sim[df_jaccard_sim['JaccardSimilarity']>0] # remove no jaccard 

df_jaccard_sim['Jaccard_Z'] = zscore(df_jaccard_sim['JaccardSimilarity'])

df_jaccard_sim.to_csv(str('/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks/DISEASES/diseases_jaccard.csv'), index=False)


############### filter significant edges and create cytoscape network

significant_edges = df_jaccard_sim[df_jaccard_sim['Jaccard_Z'] > 1.96]
print(f"Found {len(significant_edges)} significant gene-gene interactions")

cytoscape_edges = significant_edges[['Gene1', 'Gene2', 'Jaccard_Z']].copy()
cytoscape_edges.rename(columns={
    'Gene1': 'Source',
    'Gene2': 'Target',
    'Jaccard_Z': 'Weight'
}, inplace=True)


cytoscape_edges.to_csv(str('/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks/DISEASES/diseases_filtered_network_1_96.csv'), index=False)


#################################################### MGI Phenotype

url_HMD_annotation = "https://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt"
HMD_annotation = pd.read_table(url_HMD_annotation, header=None, sep="\t")

HMD_annotation = HMD_annotation.set_axis(["Human_Marker_Symbol", "Human_Entrez_Gene_ID", "Mouse_Marker_Symbol", "MGI_Marker_Accession_ID", "Highlevel_Mammalian_Phenotype_ID", "remove"], axis=1)
HMD_annotation = HMD_annotation.drop("remove", axis = 1)
HMD_annotation.head()


url_MGI_annotation = "https://www.informatics.jax.org/downloads/reports/MGI_PhenoGenoMP.rpt"
mgi_phen_annotation = pd.read_table(url_MGI_annotation, header=None)

mgi_phen_annotation = mgi_phen_annotation.set_axis(["Allelic_Composition", "AlleleSymbols",	"Genetic_Background","Mammalian_Phenotype_ID", "PubMed_IDs", "MGI_Marker_Accession_IDs"], axis=1)
mgi_phen_annotation.head()

hp_list = list(set(mgi_phen_annotation["Mammalian_Phenotype_ID"]))
hp_rows = []

for hp in hp_list:
    mgi_ids = '|'.join(mgi_phen_annotation[mgi_phen_annotation["Mammalian_Phenotype_ID"] == hp ]["MGI_Marker_Accession_IDs"]).split("|")
    for mgi_id in mgi_ids:
        hp_rows.append([hp, mgi_id])

df_hp_mgi = pd.DataFrame(hp_rows, columns=['Phenotype', 'MGI_id']) ## phenotype - gen dataframe

df_hp_mgi

##### merge the MGI id with the gene symbol. 

df_hp_mgi = df_hp_mgi.merge(HMD_annotation.loc[:,("Human_Marker_Symbol", "MGI_Marker_Accession_ID")], left_on='MGI_id', right_on='MGI_Marker_Accession_ID', how='left')

# from collections import Counter
# counts = Counter(HMD_annotation.MGI_Marker_Accession_ID)
# duplicates = [item for item, count in counts.items() if count > 1]

################################
## numero de genes únicos
unique_genes = list(set(df_hp_mgi["Human_Marker_Symbol"].dropna()))


print("Hay ", str(len(unique_genes)), " genes únicos.")

#######################################33

gene_disease_dict = {}

for gen in unique_genes:
    gene_disease_dict[gen] = list(df_hp_mgi[df_hp_mgi["Human_Marker_Symbol"] == gen]["Phenotype"])

print("En total, hay " + str(len(gene_disease_dict.keys())) + " genes únicos diferentes")
## En total, hay 14646 genes únicos diferentes

##### calculate jaccard between all pairs

similarities = {}

for gene1, gene2 in combinations(gene_disease_dict.keys(), 2):
    sim = jaccard_similarity(set(gene_disease_dict[gene1]), set(gene_disease_dict[gene2]))
    similarities[(gene1, gene2)] = sim


# Convert to a DataFrame
df_jaccard_sim = pd.DataFrame(
    [(g1, g2, sim) for (g1, g2), sim in similarities.items()],
    columns=['Gene1', 'Gene2', 'JaccardSimilarity']
)


####### calculate zscore

df_jaccard_sim = df_jaccard_sim[df_jaccard_sim['JaccardSimilarity']>0] # remove no jaccard 

df_jaccard_sim['Jaccard_Z'] = zscore(df_jaccard_sim['JaccardSimilarity'])

df_jaccard_sim.to_csv(str('/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks/MGI_phenotype/MGIPhenotype_jaccard.csv'), index=False)


############### filter significant edges and create cytoscape network

significant_edges = df_jaccard_sim[df_jaccard_sim['Jaccard_Z'] > 1.96]
print(f"Found {len(significant_edges)} significant gene-gene interactions")

cytoscape_edges = significant_edges[['Gene1', 'Gene2', 'Jaccard_Z']].copy()
cytoscape_edges.rename(columns={
    'Gene1': 'Source',
    'Gene2': 'Target',
    'Jaccard_Z': 'Weight'
}, inplace=True)


cytoscape_edges.to_csv(str('/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks/MGI_phenotype/MGIPhenotype_network_1_96.csv'), index=False)






















