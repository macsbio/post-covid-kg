#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 10:25:33 2024

@author: alejandroadriaquelozano
"""
# Set path
import os 
from pathlib import Path
path = Path(__file__).resolve().parent
os.chdir(path)

from neo4j import GraphDatabase
import pandas as pd
import networkx as nx
import py4cytoscape as p4c
import math


# Connect to Neo4j database
uri = "bolt://localhost:7687"
username = "neo4j"
password = "123456789"
driver = GraphDatabase.driver(uri, auth=(username, password))

# Define Neo4j query ADD SCORES!!!
neo4j_query = """
MATCH (g1:gene)-[r1:stringdb_link_to]-(g2:gene)
RETURN g1.name AS gene1, g2.name AS gene2, null AS gene, r1.type AS interactionType, r1.score as score ,null AS drug

UNION

MATCH (gene:gene)-[interaction:inhibits|activates]->(drug:drug)
RETURN null AS gene1, null AS gene2, gene.name AS gene, interaction.type AS interactionType, null as score ,drug.name AS drug
"""

# Execute Neo4j query
with driver.session() as session:
    result = session.run(neo4j_query)
    data = [record for record in result]
 
data_df=pd.DataFrame(data)
data_df = data_df.rename(columns={0: 'gene1', 1: 'gene2',2:'gene',3:'interaction',4:'score',5:'drug'})

data_ppi = data_df.dropna(subset=['gene2'])
#filter out relationships in STRING by score < 0.4, 0.7 and 0.9...
data_ppi = data_ppi[data_ppi['score'] >= 0.9]


data_drugs = data_df.dropna(subset=['drug'])


# P-value from the file
dea= pd.read_excel('BulkRNA-Seq-DEGs.xlsx')

dea = dea[dea['p-Value'] < 0.05]
dea = dea[['gene','p-Value','logFC']]
dea['gene'] = dea['gene'].astype(str)
dea['p-Value'] = dea['p-Value'].astype(float)

# Create a graph from DataFrame

drugs = set(data_drugs['drug'])
genes_1 = set(data_drugs['gene'])
genes_2 = set(data_ppi['gene1'])
genes_3 = set(data_ppi['gene2'])

genes= genes_1.union(genes_2, genes_3)
dea = dea[dea['gene'].isin(genes)]

G_ppi = nx.from_pandas_edgelist(data_ppi, 'gene1', 'gene2','score')
G_drugs = nx.from_pandas_edgelist(data_drugs, 'gene', 'drug')

G_edges = nx.compose(G_ppi,G_drugs)

G_nodes = nx.Graph()
G_nodes.add_nodes_from(drugs, node_type='drug')

for index, row in dea.iterrows():
    gene = row['gene']
    p_value = row['p-Value']
    logFC =  row['logFC']
    G_nodes.add_node(gene,p_value=p_value,logFC=logFC,node_type='gene')

G = nx.compose(G_nodes, G_edges)

'''
output_path='ppi_network_09.graphml'


nx.write_graphml(G, output_path, named_key_ids=True)
p4c.import_network_from_file(output_path)

# Create a visual style 
style_name = 'PPI_drus'
if style_name not in p4c.get_visual_style_names():
    p4c.create_visual_style(style_name)

# Apply the style
p4c.set_visual_style(style_name)

# Set node size mapping, shapes and colours
values=['gene','drug']
shapes=['ELLIPSE', 'TRIANGLE']
sizes= [40,50]
p4c.set_node_size_mapping('node_type',values,sizes=sizes,mapping_type='d')
p4c.set_node_shape_mapping('node_type', values, shapes)

p4c.set_node_color_default('#228B22')


control_points = [-3, 0, 3]
colors = ['#1f78b4', '#999999', '#e31a1c']

p4c.set_node_color_mapping('logFC',control_points,colors)
'''

# we are going to retrieve info of drug targets
drug_targets = data_drugs.drop(columns=['gene1','gene2','score'])
drug_targets = pd.merge(drug_targets,dea,on='gene',how='left')

degree_targets = drug_targets.groupby('gene')['drug'].agg(list).reset_index()

degree_targets['degree'] = degree_targets['drug'].apply(len)

import matplotlib.pyplot as plt

# Using the 'count' column from the previous DataFrame to plot the degree distribution
counts = degree_targets['degree']

plt.figure(figsize=(8, 6))
plt.hist(counts, bins=range(1, counts.max() + 2), align='left', rwidth=0.8)
plt.title('Degree Distribution')
plt.xlabel('Degree')
plt.ylabel('Frequency')
plt.xticks(range(1, counts.max() + 1))
plt.grid(axis='y', alpha=0.75)
plt.show()

average_degree = degree_targets['degree'].mean()
activation_counts = drug_targets['interaction'].value_counts()

def categorize_value(val):
    if val > 1:
        return 1
    elif val < -1:
        return -1
    else:
        return 0

# Apply the function to the 'numbers' column and create a new column
drug_targets_no_duplicates = drug_targets.drop_duplicates(subset=['gene'])

drug_targets_no_duplicates['activation_state'] = drug_targets_no_duplicates['logFC'].apply(categorize_value)

logFC_count=drug_targets_no_duplicates['activation_state'].value_counts()



