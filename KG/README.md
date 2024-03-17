# Post-COVID Knowledge Graph (KG) 🦠🔍📊

## Knowledge Graph Construction 🧩
In this folder you will find the data and scripts used for constructing the post-COVID-19 Knowledge graph.
First of all be sure you have installed pyBiodataFuse in your environment, otherwise install using:

```bash
   $ pip install git+https://github.com/BioDataFuse/pyBiodatafuse.git
 ```
Input files:
- input_genes.csv: this file contains the genes you want to use as source nodes
- post-COVID-19_genes.csv: this is a set of genes with a known relationship with post-COVID-19 according literature. Only use if your project has a focus on post-COVID-19

Modified Scripts: those are scripts from the package pyBiodataFuse that have been modified for this project, they enable a smooth generation of the graph
- drug_diseease_annotator.py
- new_disgenet_annotator.py
- opentargets.py
- minerva.py
- generator.py
- neo4j_exporter.py

Script: this is the file that generates the KG
- gene_workflow_generator.py

Output files:
- graph.csv: formatted table with all the information, however not in triplets format
- graph.graphml: graph ready for importing into neo4j
- combined_df.plk: graph in plk format
