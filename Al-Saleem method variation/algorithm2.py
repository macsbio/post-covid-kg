
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

# Drugs experimental OR g.logFC < 1 positive and downregulated genes 
# TODO add also lower? What to do with inf values?
neo4j_query = """
MATCH (g:gene)-[r:part_of_go]->(bp:`gene ontology` {go_type: 'P'})
WHERE g.logFC > 1 OR g.logFC < -1
WITH bp, collect(DISTINCT g) AS RegulatedGenes
WITH bp, RegulatedGenes, size(RegulatedGenes) AS RegulatedGenesCount
WHERE RegulatedGenesCount >= 4
UNWIND RegulatedGenes AS Genes
MATCH (Genes)-[:activates|inhibits]->(d:drug)
WITH d, collect(DISTINCT Genes) AS AssociatedGenes, collect(DISTINCT bp.name) AS Ontologies
RETURN d.name, d.adverse_effect_count as SideEffect, [g IN AssociatedGenes | g.name] AS GeneNames, Ontologies
"""


with driver.session() as session:
    result = session.run(neo4j_query)
    data = [record for record in result]
 
drugs_exp=pd.DataFrame(data)
drugs_exp=drugs_exp.rename(columns={1:'side_effects',2:'genes_exp',3:'bp_exp'})

#Literature drugs
neo4j_query = """
// Find genes associated with Post-COVID-19 and count their associated biological processes
MATCH (g:gene)-[:associated_with]->(:disease{name:'Post-COVID-19'})
WITH g
MATCH (g)-[:part_of_go]->(bp:`gene ontology`{go_type:'P'})
WITH bp, COUNT(g) AS GeneCount
WHERE GeneCount >= 4
WITH collect(bp.name) AS RelevantBPs

// Find genes related to those relevant BP 
MATCH (g:gene)-[:part_of_go]->(bp:`gene ontology`{go_type:'P'})
WHERE bp.name IN RelevantBPs
MATCH (g)-[:activates|inhibits]->(d:drug)
WITH RelevantBPs, d, collect(g) AS Genes

// Return drugs, their side effects, associated genes, and biological processes present in 4 or more genes
RETURN d.name AS DrugName, [gene IN Genes | gene.name] as Genes, RelevantBPs
"""

with driver.session() as session:
    result = session.run(neo4j_query)
    data = [record for record in result]
 
drugs_lit=pd.DataFrame(data)
drugs_lit=drugs_lit.rename(columns={1:'genes_lit',2:'bp_lit'})



#Combine dataframes (ex and literature)
drugs_df= pd.merge(drugs_exp, drugs_lit,on=0,how='outer')

#### Score calculation

#GR
neo4j_query = """
MATCH (g:gene)
WITH collect(g.name) AS geneNames
UNWIND geneNames AS geneName
MATCH (g:gene {name: geneName})-[r:activates|inhibits]->(d:drug)
WITH geneName, apoc.coll.toSet(collect(d)) AS uniqueDrugs
WITH geneName, size(uniqueDrugs) AS drugCount, toFloat(1) / log(1.5+size(uniqueDrugs)) AS geneSig
RETURN geneName, geneSig
"""
    
with driver.session() as session:
      result = session.run(neo4j_query)
      data = [record for record in result]
gene_rarity=pd.DataFrame(data)  

gene_rarity=gene_rarity.rename(columns={1:'GR'})

#Biological process rarity (BPR)
neo4j_query="""
MATCH (b:`gene ontology`) 
WHERE b.go_type = 'P'
WITH collect(b) AS bs
UNWIND bs AS b
MATCH (b)<-[:part_of_go]-(g:gene)-[r:activates|inhibits]->(d:drug)
WITH b, size(apoc.coll.toSet(COLLECT(d))) AS countUniqueDrugs, toFloat(1) / log(1.5+size(apoc.coll.toSet(COLLECT(d)))) AS bpSig
RETURN b.name, bpSig
"""
with driver.session() as session:
      result = session.run(neo4j_query)
      data = [record for record in result]
bp_rarity=pd.DataFrame(data)  

bp_rarity=bp_rarity.rename(columns={1:'BPR'})

#Cytokine release
neo4j_query="""
MATCH (d:drug)<-[:activates|inhibits]-(g:gene)-[:part_of_go]->(:`gene ontology`{name:'positive regulation of neuroinflammatory response'})
return d.name
"""
with driver.session() as session:
      result = session.run(neo4j_query)
      data = [record for record in result]
neuroinflammatory =pd.DataFrame(data)  

neuroinflammatory['neuroinflammatory response']=1



# Fix drugs
drugs_df = pd.merge(drugs_df,neuroinflammatory,on=0,how='left')
drugs_df['neuroinflammatory response'] = drugs_df['neuroinflammatory response'].fillna(0)
drugs_df['side_effects'] = drugs_df['side_effects'].fillna(0)

drugs_df['side_effects'] = drugs_df['side_effects'] +1.5 #TODO Does this makes sense?



#Score function

results= pd.DataFrame()

for index,row in drugs_df.iterrows():
    side_effects= math.log(row['side_effects'])
    gene= row[0]
    
    sum_gr_exp=0
    sum_gr_lit=0
    sum_bpr_exp=0
    sum_bpr_lit=0
   
    
    if type(row['genes_exp'])!= float:
        for gene in row['genes_exp']:
            gr_exp= gene_rarity[gene_rarity[0] == gene]
            gr_exp = gr_exp['GR'].iloc[0]
            sum_gr_exp=gr_exp+sum_gr_exp
        
        
    if type(row['genes_lit'])!= float:
        for gene in row['genes_lit']:
            gr_lit= gene_rarity[gene_rarity[0] == gene]
            gr_lit = gr_lit['GR'].iloc[0]
            sum_gr_lit=gr_lit+sum_gr_lit
        
            
    if type(row['bp_exp'])!= float:
       for bp in row['bp_exp']:
           bpr_exp= bp_rarity[bp_rarity[0] == bp]
           bpr_exp = bpr_exp['BPR'].iloc[0]
           sum_bpr_exp=bpr_exp+sum_bpr_exp


    if type(row['bp_lit'])!= float:
       for bp in row['bp_lit']:
           bpr_lit= bp_rarity[bp_rarity[0] == bp]
           bpr_lit = bpr_lit['BPR'].iloc[0]
           sum_bpr_lit=bpr_lit+sum_bpr_lit

    
    numerator = (sum_gr_lit+sum_bpr_lit)+(sum_gr_exp+sum_bpr_exp)
    score= (numerator/side_effects)+2*row['neuroinflammatory response']
    dict_entry={'drug':row[0],'score':score}
    df_dictionary = pd.DataFrame([dict_entry])
    results = pd.concat([results, df_dictionary], ignore_index=True)
    print(score, '\t',row[0])
    
'''
results=results.sort_values(by='score', ascending=False)
results2= results.head(30)
results2.to_csv('/Users/alejandroadriaquelozano/Documents/Internships/MacsBio/results2.csv')


##Common values between results
results104= pd.read_csv('/Users/alejandroadriaquelozano/Documents/Internships/MacsBio/results1_04.csv')
results107= pd.read_csv('/Users/alejandroadriaquelozano/Documents/Internships/MacsBio/results1_07.csv')
results109= pd.read_csv('/Users/alejandroadriaquelozano/Documents/Internships/MacsBio/results1_09.csv')

results104=results104.head(30)
results107=results107.head(30)
results109=results109.head(30)

common_values04 = list(set(results104['name']) & set(results2['drug']))
common_values07 = list(set(results107['name']) & set(results2['drug']))
common_values09 = list(set(results109['name']) & set(results2['drug']))
'''

#Count gene list for figure

genes_lit=list()
for x in drugs_df['genes_lit']:
    if type(x) == list:
        for y in x: 
            if type(y)==str:
                genes_lit.append(y)
genes_lit=set(genes_lit)



genes_exp=list()
for x in drugs_df['genes_exp']:
    if type(x) == list:
        for y in x: 
            if type(y)==str:
                genes_exp.append(y)
genes_exp=set(genes_exp)

intersected_set = genes_exp.intersection(genes_lit)
intersected_set2 = genes_lit.intersection(genes_exp)

