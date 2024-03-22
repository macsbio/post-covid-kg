# Post-COVID Knowledge Graph (KG) 🦠🔍📊

## Data Exploration and Visualization 📊👀

### Network visalisation in cytoscape
In this part we create a cytoscape network of the gene-gene and gene-drug interctions using the graph in neo4j. 
First, be sure you have your graph in neo4j.

Script 
- ppy_cytoscape.py: this criptos connects with your neo4j graph and queries it to extract  gene-gene and gene-drug information. After that creates graphml networks and imports it to cytoscape.

Output
- 3 grapgml networksusing 3 different STRING scores for gene-gene interactions
- ppi_results.cys: cytoscape session where you can visualise the networks
