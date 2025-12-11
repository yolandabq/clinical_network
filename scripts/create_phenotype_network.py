#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 15:56:24 2025

@author: yolanda

Script for merging HPO, MGI and Diseases network to create the clinical network, and the Intellectual Disability subnetwork using the genes from PanelApp (Green and Ambar).

"""

import networkx as nx
import pandas as pd
import numpy as np

######### create functions 

def get_network_info(G): 
    degree_dict = dict(G.degree())
    top_nodes = sorted(degree_dict.items(), key=lambda x: x[1], reverse=True)[:10]
    print("Top connected nodes:", top_nodes)
    print("Density:", nx.density(G))
    print("Nodes:", G.number_of_nodes())
    print("Edges:", G.number_of_edges())

    # Número de componentes
    num_componentes = nx.number_connected_components(G)

    # Lista de componentes (cada uno es un conjunto de nodos)
    componentes = list(nx.connected_components(G))

    # Tamaños de cada componente
    sizes = [len(c) for c in componentes]

    print("Number of connected components:", num_componentes)
    print("Components size:", sizes)
    

##########

## path to individual phenotype networks:
    
files_phenotype = ["/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks/new_updated_networks/phenotype/phenotype_HPO.txt", 
         #"/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks/new_updated_networks/phenotype/phenotype_HPO_ext.txt",
         "/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks/new_updated_networks/phenotype/phenotype_MGI.txt",
         #"/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks/new_updated_networks/phenotype/phenotype_PanelApp.txt",
         "/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks/new_updated_networks/phenotype/phenotype_DISEASES.txt"]

##########
## Intellectual disability genes form PanelApp

ID_genes_file = "/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks/PanelApp/GA_panels/Intellectual_disability_GA.csv" # Green Ambar
pd_ID_genes = pd.read_table(ID_genes_file)
pd_ID_genes.head()
ID_genes = set(pd_ID_genes["gene_data.gene_symbol"])
print("Number of Intellectual Disability genes (Green-Ambar): " + str(len(ID_genes)))
#########
## Merge the three individual networks

dfs_phenotype = [pd.read_table(f, header=None, sep=" ", usecols=[0,1,3]) for f in files_phenotype]
unique_genes = [pd.unique(dfs_phenotype[x][[0,1]].values.ravel()) for x in [0,1,2]]
print("Number of unique genes in each individual network: " + str([len(unique_genes[x]) for x in [0,1,2]]))


dfs_phenotype[0][3] = dfs_phenotype[0][3].str.replace('}', '').astype(float)
dfs_phenotype[1][3] = dfs_phenotype[1][3].str.replace('}', '').astype(float)
dfs_phenotype[2][3] = dfs_phenotype[2][3].str.replace('}', '').astype(float)

# Step 1: concat the three dataframes
df_all = pd.concat(dfs_phenotype, ignore_index=True)

# Step 2: sort the pair of genes son (A,B) and (B,A) count as the same interaction
df_all['gene_pair'] = df_all.apply(
    lambda row: tuple(sorted([row[0], row[1]])), axis=1
)

# Step3: We keep the max weight for each interacction
df_max = df_all.groupby('gene_pair')[3].max().reset_index()

# Step 4: We create columns Source and Target
df_max[['Source', 'Target']] = pd.DataFrame(df_max['gene_pair'].tolist(), index=df_max.index)

# Step 4: We only keep Source, Target and Weight columns
df_final_clinical_nx = df_max[['Source', 'Target', 3]].rename(columns={3: 'weight'})
df_final_clinical_nx.head()

df_final_clinical_nx.to_csv(str('/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks/new_updated_networks/phenotype/phenotype_merged_networks.csv'), index=False)


# Step 5: We create the networkx object
    
G = nx.from_pandas_edgelist(df_final_clinical_nx, 'Source', 'Target', edge_attr='weight')
get_network_info(G)

## create ID subnetwork

subnetwork_nodes = [n for n in G.nodes if n in ID_genes]
print("Number of ID genes present in the clinical network: " + str(len(subnetwork_nodes)))
## G.remove_edges_from(nx.selfloop_edges(G)) # for removing self loops (A = A)
## [x for x in G.edges() if x[0] == x[1]] # there isn't self lopps

G_phenotype_ID = G.subgraph(subnetwork_nodes)

get_network_info(G_phenotype_ID)

# save the network for cytoscape
#nx.write_graphml(G_phenotype_ID, "/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks/new_updated_networks/phenotype/phenotype_network.graphml")





































