#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 17:08:41 2025

@author: yolanda
"""

import networkx as nx
import pandas as pd
import numpy as np
import igraph as ig
import leidenalg as la
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
    

def nx_to_igraph(G):
    # Nodos en orden estable
    nodes = list(G.nodes())
    node_index = {n: i for i, n in enumerate(nodes)}

    edges = [(node_index[u], node_index[v]) for u, v in G.edges()]
    weights = [G[u][v].get("weight", 1.0) for u, v in G.edges()]

    g = ig.Graph(edges=edges, directed=False)
    g.vs["name"] = nodes
    g.es["weight"] = weights

    return g, nodes


def run_leiden(g, resolution=1.0, seed=123):
    partition = la.find_partition(
        g,
        la.RBConfigurationVertexPartition,
        resolution_parameter=resolution,
        weights=g.es["weight"], 
        seed=seed
    )
    
    clusters = partition.membership
    modularity = partition.modularity
    return clusters, modularity

##########
## path to the clinical network (phenotype network)
phenotype_nx = "/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks/new_updated_networks/phenotype/phenotype_merged_networks.csv"
pd_phenotype_nx = pd.read_table(phenotype_nx, sep=",")
##########
## Intellectual disability genes form PanelApp

ID_genes_file = "/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks/PanelApp/GA_panels/Intellectual_disability_GA.csv" # Green Ambar
pd_ID_genes = pd.read_table(ID_genes_file)
pd_ID_genes.head()
ID_genes = set(pd_ID_genes["gene_data.gene_symbol"])
print("Number of Intellectual Disability genes (Green-Ambar): " + str(len(ID_genes)))
#########

# Step 1: We create the networkx object
    
G = nx.from_pandas_edgelist(pd_phenotype_nx, 'Source', 'Target', edge_attr='weight')
get_network_info(G)

## create ID subnetwork

subnetwork_nodes = [n for n in G.nodes if n in ID_genes]
print("Number of ID genes present in the clinical network: " + str(len(subnetwork_nodes)))
## G.remove_edges_from(nx.selfloop_edges(G)) # for removing self loops (A = A)
## [x for x in G.edges() if x[0] == x[1]] # there isn't self lopps

G_phenotype_ID = G.subgraph(subnetwork_nodes).copy()

get_network_info(G_phenotype_ID)

####
# isolated genes: 
list(nx.isolates(G_phenotype_ID))
# we remove the isolated genes
isolados = list(nx.isolates(G_phenotype_ID))
G_phenotype_ID.remove_nodes_from(isolados)


# Step 2: Leiden Clustering

######################
##### LEIDEN ######
####################

# we create the grapth
g_leiden, node_order = nx_to_igraph(G_phenotype_ID)

## optimization of Leiden to get best resolution value #

resolutions = np.arange(0.2, 2.01, 0.1)
results = []

for r in resolutions:
    labels, mod = run_leiden(g_leiden, resolution=r)
    n_clusters = len(set(labels))
    results.append((r, mod, n_clusters))

df_leiden = pd.DataFrame(results, columns=["resolution", "modularity", "clusters"])
df_leiden.loc[df_leiden["modularity"].idxmax()]


#### To obtain the best partition
    
best_mod = -1
best_labels = None

for s in range(200):
    labels, mod = run_leiden(g_leiden, resolution=df_leiden.loc[df_leiden["modularity"].idxmax()]["resolution"], seed=s)
    if mod > best_mod:
        best_mod = mod
        best_labels = labels

print("Mejor modularidad:", best_mod)

df_leiden_opt = pd.DataFrame({
    "gen": node_order,
    "cluster_leiden": best_labels
})


df_leiden_opt["cluster_leiden"].value_counts()

# save the clustering
df_leiden_opt.to_csv("/home/yolanda/tblab/yolanda/GLOWgenes/berta/networks_exploration/networks_exploration/louvain_clustering/leiden_clustering_2.csv", index=False, header=True)




































