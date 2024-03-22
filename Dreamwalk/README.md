# Post-COVID Knowledge Graph (KG) 🦠🔍📊

## Drug Repurposing Algorithm: DreamWalk 💊🔄
This part involves running dreamwalk algorithm in the previously generated KG. 

### Preprocessing
First, we need to preprocess the data to feed the algorithm
Input files:
- data.csv: exported data from neo4j. You just have to import the graph.graphml file and export it back in csv format.
- atc_hierarchy: atc hierarchy for drugs ibtained using the drug_parser.py in the utils folder.
- DREAMwalk folder: contaings the algorithm scripts

Scripts:
- preprocess_dreamwalk.py

Output file:
- preprocessed_graph.csv: this is going to be used as input file for the algorithm

### Algorithm implementation
Input files:
-  preprocessed_graph.csv:

Scripts:
- dreamwalk_script.py

Output files:
- graph.txt
- nodetypes.tsv
- dis_sim.tsv
- similarity_graph_drugs.tsv
- similarity_graph.txt
- hierarchy.csv
- embedding_file.plk
- dda_files folder: drug-disease interaction files randomly generted
- results folder: contains the models trained with the data

### Alternative results
- 0.9 Workflow: contains the same workflow but using 0.9 as cutoff for gene-gene interactions baed on STRING scores. 
