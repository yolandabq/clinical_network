#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 11:14:08 2024

@author: yolanda

Code for obtaining panelAPP genes using the API
"""

import requests
import json
import pandas as pd

out_dir = "/mnt/tblab/yolanda/GLOWgenes/berta/networks/PanelApp/"

# primero, obtenemos la lista de paneles
server = "https://panelapp.genomicsengland.co.uk"
#ext = f"/api/v1/panels/signedoff/"
ext = f"/api/v1/panels/"

r = requests.get(server + ext, headers={"Content-Type": "application/json"})

expected_panels = r.json()["count"]
# df columns: 'Name', 'DiseaseSubGroup', 'DiseaseGroup', 'CurrentVe+rsion',
# 'CurrentCreated', 'Number_of_Genes', 'Number_of_STRs', 'Number_of_Regions',
# 'Panel_Id', 'Relevant_disorders', 'Status', 'PanelTypes'

GEL_panel_app_df = pd.json_normalize(r.json(), record_path=["results"])

# Reiterate over remaining pages of data
while  r.json()["next"] is  not  None:
	r = requests.get(
	r.json()["next"], headers={"Content-Type": "application/json"})
	GEL_panel_app_df = pd.concat([GEL_panel_app_df, pd.json_normalize(r.json(), record_path=["results"])], ignore_index=True)


ids_panel=(GEL_panel_app_df[["id", "name", "version", "disease_group", "types"]])
#ids_panel.to_csv('/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/panels_class.tsv', sep="\t", index=False)

# ahora obtenemos cada panel
server = "https://panelapp.genomicsengland.co.uk"
#ext = f"/api/v1/panels/{panel_id}/genes/?version={panel_version}"


for ind in ids_panel.index:

    
    panel_id = ids_panel['id'][ind]
    panel_version = ids_panel['version'][ind]
    panel_name = ids_panel['name'][ind]
    panel_name = panel_name.replace("/", "_")
    
    ext = f"/api/v1/panels/{panel_id}/genes/"
    #print(f"{server}{ext}")
    #print(panel_name)
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})
    
    expected_genes = r.json()["count"]
    #print(expected_genes)

    GEL_panel_data = pd.json_normalize(r.json(), record_path=["results"])
    
    while  r.json()["next"] is  not  None:
    	r = requests.get(
    	r.json()["next"], headers={"Content-Type": "application/json"})
    	GEL_panel_data = pd.concat([GEL_panel_data, pd.json_normalize(r.json(), record_path=["results"])], ignore_index=True)
        
    #print(GEL_panel_data.shape[0])
    
    if expected_genes != GEL_panel_data.shape[0]:
        print("revisar panel!!!!!!!!!!! " + panel_name)
        
    if expected_genes == 0: 
        print("REVISAR PANEL CON 0 GENES: " + panel_name)
        
    else: 
        GEL_panel_data=GEL_panel_data.rename(columns={"gene_data.hgnc_symbol": "gene"})
        GA_panel = GEL_panel_data[(GEL_panel_data["confidence_level"] == "3") | (GEL_panel_data["confidence_level"] == "2") ]
        #"_".join(panel_name.split(" "))
        GEL_panel_data.to_csv(str(out_dir + "_".join(panel_name.split(" ")) + "_GAR.csv"), index=False, sep="\t" )
        GA_panel.to_csv(str(out_dir + "_".join(panel_name.split(" ")) + "_GA.csv") , index=False, sep="\t")
