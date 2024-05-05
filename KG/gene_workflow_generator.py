#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 16:30:05 2024

@author: alejandroadriaquelozano
"""
# Set path
import os 
from pathlib import Path
path = Path(__file__).resolve().parent
os.chdir(path)

# Import modules
import pickle
import pandas as pd
from pyBiodatafuse import id_mapper
from pyBiodatafuse.annotators import stringdb, wikipathways
from pyBiodatafuse.utils import combine_sources
import generator,neo4j_exporter, new_disgenet_annotator,drug_disease_annotator,minerva,opentargets

# Differentially expressed genes in post-COVID-19
data_input = pd.read_csv('input_genes.csv')


data_input.head()

bridgdb_df, bridgdb_metadata = id_mapper.bridgedb_xref(
    identifiers=data_input,
    input_species="Human",
    input_datasource="HGNC",
    output_datasource="All",
)
bridgdb_df.head()

#Diseases Disgenet
disgenet_result = new_disgenet_annotator.get_disgenet_diseases(bridgdb_df)

#Location
loc_df, opentargets_loc_metadata = opentargets.get_gene_location(bridgedb_df=bridgdb_df)
loc_df.head()

#GO
go_process_df, opentargets_go_metadata = opentargets.get_gene_go_process(bridgedb_df=bridgdb_df)
go_process_df.head()

#reactome pathways
reactome_process_df, opentargets_process_metadata = opentargets.get_gene_reactome_pathways(
    bridgedb_df=bridgdb_df
)
reactome_process_df.head()

# Drugs
drug_df, opentargets_drug_metadata = opentargets.get_gene_drug_interactions(bridgedb_df=bridgdb_df)
drug_df.head()

'''Not interested in these
Diseases (opentargets)
disease_df, opentargets_disease_metadata = opentargets.get_gene_disease_associations(
    bridgedb_df=bridgdb_df
)
disease_df.head()
'''

#Wikipathways
wp_df, wp_metadata = wikipathways.get_gene_wikipathways(bridgedb_df=bridgdb_df)
wp_df.head()

#Protein-protein STRING
ppi_df, ppi_metadata = stringdb.get_ppi(bridgedb_df=bridgdb_df)
ppi_df.head()

#MINERVA
project_list_df = minerva.list_projects()
map_components = minerva.get_minerva_components(project_list_df, map_name='COVID19 Disease Map', get_reactions=False )
minerva_df = minerva.get_gene_minerva_pathways(bridgdb_df,map_components)

# DRUG DISEASE RELATIONSHIPS
drug_disease = drug_disease_annotator.get_drug_disease_interactions(drug_df,disgenet_result)

# Adding post-COVID-19 disease nodes using csv file (only run this code if you are interested in post-COVID-19)
post_covid_genes = pd.read_csv('post-COVID-19_genes.csv')
row=0
post_covid_node =  {'disease_name': 'Post-COVID-19', 'diseaseid': 'C0000000'}
for x in disgenet_result['identifier']:
    if x.strip() in list(post_covid_genes['identifier']):
        lol=(disgenet_result['DisGeNET'][row])
        disgenet_result['DisGeNET'][row].append(post_covid_node)
        
    row=row+1
        

combined_df = combine_sources(
    [
        disgenet_result,
        loc_df,
        go_process_df,
        reactome_process_df,
        drug_df,
        #disease_df,
        wp_df,
        ppi_df,
        minerva_df
    ]
)

#NCBI IDs (only if you are interested in also obtaining NCBI IDs)
ncbi_genes = bridgdb_df[bridgdb_df['target.source'] == 'NCBI Gene']
ncbi_genes = ncbi_genes.drop(columns=['target.source','identifier.source'])
ncbi_genes = ncbi_genes.rename(columns={'target': 'ncbi_id'})

combined_df= pd.merge(combined_df, ncbi_genes,on='identifier',how='left')


#Export formatted table
combined_df.to_csv('graph.csv')


#Export graph
with open("combined_df.pkl", "wb") as out:
    pickle.dump(combined_df, out)
    
combined_df = generator.load_dataframe_from_pickle("combined_df.pkl")

pygraph = generator.generate_networkx_graph(combined_df,drug_disease)

neo4j_exporter.save_graph_to_neo4j_graphml(pygraph, 'graph.graphml')



