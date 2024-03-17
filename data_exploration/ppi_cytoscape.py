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


# Connect to Neo4j database
uri = "bolt://localhost:7687"
username = "neo4j"
password = "123456789"
driver = GraphDatabase.driver(uri, auth=(username, password))

# Define Neo4j query ADD SCORES!!!
neo4j_query = """
MATCH (g1:gene)-[r1:interacts_with]-(g2:gene)
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
#filter out relationships in STRING by score < 0.7 and 0.9...
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
    G_nodes.add_node(gene,p_value=p_value,logFC=logFC)

G = nx.compose(G_nodes, G_edges)


output_path='ppi_network_09.graphml'


nx.write_graphml(G, output_path, named_key_ids=True)
p4c.import_network_from_file(output_path)


#Applpy style 
values=['gene','drug']
shapes=['ELLIPSE', 'TRIANGLE']
p4c.set_node_shape_mapping('node_type', values, shapes)


column = 'p-value'
control_points = [0.0, 0.05]
colors = ['#FFFFFF', '#DD8855']

p4c.set_node_color_mapping('p_value',control_points,colors)
    

