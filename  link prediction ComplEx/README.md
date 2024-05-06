# Post-COVID Knowledge Graph (KG) 🦠🔍📊

## Drug Repurposing Algorithm: Trouillon method variation (ComplEx) 💊🔄
This part involves running ComplEx algorithm in the previously generated post-COVID-19 KG. (https://proceedings.mlr.press/v48/trouillon16.html)

### Algorithm implementation
There's no need for an input archive because we are querying directly from Neo4j. Just make sure your KG is correctly implemented in neo4j.

Scripts:
- link_prediction.py

Output files:
- ComplEx_results.csv
- filtered_ranks and raw_ranks for the testing set -> works as quality measure of the graph predictions
- Accuracy graph of the model train set vs the validation set
