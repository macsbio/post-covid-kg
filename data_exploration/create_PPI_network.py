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
RETURN g1.name AS gene1, g2.name AS gene2, r1.name AS interactionType, r1.score as score
"""

# Execute Neo4j query
with driver.session() as session:
    result = session.run(neo4j_query)
    data = [record for record in result]
 
data_df=pd.DataFrame(data)
data_df = data_df.rename(columns={0: 'gene1', 1: 'gene2',2:'interaction',3:'score'})

data_ppi = data_df.dropna(subset=['gene2'])
#filter out relationships in STRING by score < 0.4, 0.7 and 0.9...
data_ppi = data_ppi[data_ppi['score'] >= 0.4]





G_ppi = nx.from_pandas_edgelist(data_ppi, 'gene1', 'gene2','score')

genes1 =set(data_df['gene1'])
genes2= set(data_df['gene2'])

genes = list(genes1.intersection(genes2))

G_nodes = nx.Graph()
G_nodes.add_nodes_from(genes, node_type='gene')


G = nx.compose(G_nodes, G_ppi)


output_path='ppi_network_04.graphml'


nx.write_graphml(G, output_path, named_key_ids=True)
p4c.import_network_from_file(output_path)
