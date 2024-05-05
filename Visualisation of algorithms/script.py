#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 14:58:20 2024

@author: alejandroadriaquelozano
"""

import os 
from neo4j import GraphDatabase
from pathlib import Path
from pyvis.network import Network


path = Path(__file__).resolve().parent
os.chdir(path)

import pandas as pd

results2 = pd.read_csv('results_3.csv')


# Connect to Neo4j database
uri = "bolt://localhost:7687"
username = "neo4j"
password = "123456789"
driver = GraphDatabase.driver(uri, auth=(username, password))



drug_list = results2['name'].head(10).tolist()


query = """
MATCH (d:drug)<-[:activates|inhibits]-(g:gene)-[:part_of_pathway]->(p:pathway)
WHERE d.name IN $drugs
RETURN d.name AS Drug, g.name AS Gene, p.name AS Pathway
"""

with driver.session() as session:
    result = session.run(query, drugs=drug_list)
    data = [record for record in result]


# Visualisation

net = Network(notebook=True)

# Add nodes and edges
for record in data:
    net.add_node(record['Drug'], label=record['Drug'], color='red', size=20)
    net.add_node(record['Gene'], label=record['Gene'], color='green', size=15)
    net.add_node(record['Pathway'], label=record['Pathway'], color='blue', size=10)
    net.add_edge(record['Drug'], record['Gene'], title='AFFECTS')
    net.add_edge(record['Gene'], record['Pathway'], title='PART_OF')

# Generate network
net.show("example.html")


