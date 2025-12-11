#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 29 12:05:31 2025

@author: yolanda
"""

import networkx
import pandas as pd
import obonet
from itertools import combinations
from scipy.stats import zscore


#### functions

def jaccard_similarity(set1, set2):
    # Efficient Jaccard similarity calculation
    intersection = len(set1 & set2)
    union = len(set1 | set2)
    return intersection / union if union != 0 else 0

######### we select hpo terms from Phenotypic abnormality HP:0000118

#https://pypi.org/project/obonet/

# Read the hpo ontology
url = 'https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2025-05-06/hp-base.obo'
graph = obonet.read_obo(url)

# Number of nodes
print("The HPO ontology has "+ str(len(graph)) + " terms") 

# Check if the ontology is a DAG
networkx.is_directed_acyclic_graph(graph)

phen_abnormality_hpos = list(networkx.ancestors(graph, "HP:0000118"))
len(phen_abnormality_hpos) # children of phenotypic abnormality

######################## 

## read the annotation file

annotation_file_options = ["genes_to_phenotype", "phenotype_to_genes"] # or phenotype_to_genes

for annotation_file in annotation_file_options:

    url_hpo_annotation = str('https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2025-05-06/' + annotation_file + '.txt')
    hpo_annotation = pd.read_table(url_hpo_annotation)
    
    # filter only hpo under phenotypic abnormality
    
    hpo_annotation = hpo_annotation[hpo_annotation.hpo_id.isin(phen_abnormality_hpos)]
    hpo_annotation = hpo_annotation[hpo_annotation["gene_symbol"] != "-"]
    ##### create gen2phen_dict
    
    gen2phen_dict = {}
    
    for gen in list(set(hpo_annotation.gene_symbol)):
        gen2phen_dict[gen] = list(hpo_annotation[hpo_annotation["gene_symbol"] == gen]["hpo_id"])
    
    print("En total, hay " + str(len(gen2phen_dict.keys())) + " genes únicos diferentes")
    ## En total, hay 5124 genes únicos diferentes
    
    ##### calculate jaccard between all pairs
    
    similarities = {}
    
    for gene1, gene2 in combinations(gen2phen_dict.keys(), 2):
        sim = jaccard_similarity(set(gen2phen_dict[gene1]), set(gen2phen_dict[gene2]))
        similarities[(gene1, gene2)] = sim
    
    
    # Convert to a DataFrame
    df_jaccard_sim = pd.DataFrame(
        [(g1, g2, sim) for (g1, g2), sim in similarities.items()],
        columns=['Gene1', 'Gene2', 'JaccardSimilarity']
    )
    
    ####### calculate zscore
    
    df_jaccard_sim = df_jaccard_sim[df_jaccard_sim['JaccardSimilarity']>0] # remove no jaccard 
    
    df_jaccard_sim['Jaccard_Z'] = zscore(df_jaccard_sim['JaccardSimilarity'])
    
    df_jaccard_sim.to_csv(str('/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks/hpo/' + annotation_file +  'jaccard.csv'), index=False)
    
    
    ############### filter significant edges and create cytoscape network
    
    significant_edges = df_jaccard_sim[df_jaccard_sim['Jaccard_Z'] > 1.96]
    print(f"Found {len(significant_edges)} significant gene-gene interactions")
    
    cytoscape_edges = significant_edges[['Gene1', 'Gene2', 'Jaccard_Z']].copy()
    cytoscape_edges.rename(columns={
        'Gene1': 'Source',
        'Gene2': 'Target',
        'Jaccard_Z': 'Weight'
    }, inplace=True)
    
    
    cytoscape_edges.to_csv(str('/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks/hpo/' + annotation_file +  '_1_96.csv'), index=False)



























































