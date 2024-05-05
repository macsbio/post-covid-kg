#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 12:25:01 2024

@author: alejandroadriaquelozano
"""
import matplotlib.pyplot as plt
from math import isclose
from sklearn.decomposition import PCA
import os
import networkx as nx
import numpy as np
import pandas as pd
from stellargraph import StellarGraph, datasets
from stellargraph.data import EdgeSplitter
from collections import Counter
import multiprocessing
from IPython.display import display, HTML
from sklearn.model_selection import train_test_split
from stellargraph.data import BiasedRandomWalk
from gensim.models import Word2Vec
from stellargraph import datasets, utils
from tensorflow.keras import callbacks, optimizers, losses, metrics, regularizers, Model
import numpy as np
import pandas as pd

from stellargraph.mapper import KGTripleGenerator
from stellargraph.layer import ComplEx

from IPython.display import HTML

import os 
import multiprocessing
from pathlib import Path
path = Path(__file__).resolve().parent
os.chdir(path)


# 
data_all= pd.read_csv('data.csv')

nodes= data_all.dropna(subset=['_id'])
relationships = data_all.dropna(subset=['_start'],ignore_index=True)

drugs=nodes[nodes['_labels']==':drug']
drugs_df = drugs[['_id']].set_index('_id')


diseases=nodes[nodes['_labels']==':disease']
diseases_df= diseases[['_id']].set_index('_id')

genes=nodes[nodes['_labels']==':gene']
genes_df= genes[['_id','logFC']].set_index('_id')



types_of_interest = ['stringdb_link_to', 'activates', 'inhibits', 'associated_with', 'treated_with']

# Filter the DataFrame based on whether the '_type' column's values are in the above list
relationships = relationships[relationships['_type'].isin(types_of_interest)]


relationships_df_2= relationships[['_start', '_end','_type']].reset_index(drop=True)
relationships_df_2= relationships_df_2.rename(columns={'_start':'source','_end':'target','_type':'label'})
nerworder=['source','label','target']
relationships_df_2=relationships_df_2[nerworder].reset_index(drop=True)

graph = StellarGraph({"drug": drugs_df, "disease": diseases_df,"gene":genes_df}, relationships_df_2,edge_type_column='label')

print(graph.info())



#train test
train,test_valid= train_test_split(relationships_df_2,test_size=0.1)
test,valid =train_test_split(test_valid,test_size=0.5)


#ComplEx method
epochs = 50
embedding_dimension = 200
negative_samples = 10

gen = KGTripleGenerator(
    graph, batch_size=len(train) // 100  # ~100 batches per epoch
)

wn18_complex = ComplEx(
    gen,
    embedding_dimension=embedding_dimension,
    embeddings_regularizer=regularizers.l2(1e-7),
)

inp, out = wn18_complex.in_out_tensors()

model = Model(inputs=inp, outputs=out)

model.compile(
    optimizer=optimizers.Adam(learning_rate=0.001),
    loss=losses.BinaryCrossentropy(from_logits=True),
    metrics=[metrics.BinaryAccuracy(threshold=0.0)],
)

train_gen = gen.flow(
    train, negative_samples=negative_samples, shuffle=True
)
valid_gen = gen.flow(valid, negative_samples=negative_samples)

es = callbacks.EarlyStopping(monitor="val_loss", patience=10)
history = model.fit(
    train_gen, validation_data=valid_gen, epochs=epochs, callbacks=[es]
)

utils.plot_history(history)

raw_ranks_test, filtered_ranks_test = wn18_complex.rank_edges_against_all_nodes(
    gen.flow(test), graph
)

# helper function to compute metrics from a dictionary of name -> array of ranks
def results_as_dataframe(name_to_results):
    return pd.DataFrame(
        name_to_results.values(),
        columns=["mrr", "hits at 1", "hits at 3", "hits at 10"],
        index=name_to_results.keys(),
    )


def summarise(name_to_ranks):
    return results_as_dataframe(
        {
            name: (
                np.mean(1 / ranks),
                np.mean(ranks <= 1),
                np.mean(ranks < 3),
                np.mean(ranks <= 10),
            )
            for name, ranks in name_to_ranks.items()
        }
    )

## Predict post covid drugs

disease_drug_pairs = pd.DataFrame({
    'source': [18920] * len(drugs_df),
    'label': ['treated_with'] * len(drugs_df),
    'target': drugs_df.index.tolist()
})

raw_ranks, filtered_ranks = wn18_complex.rank_edges_against_all_nodes(
    gen.flow(disease_drug_pairs), graph
)

## Process datasets
raw_ranks_df= pd.DataFrame(raw_ranks)
filtered_ranks_df=pd.DataFrame(filtered_ranks)
results=pd.DataFrame()
results['_id']=disease_drug_pairs['target']
results['raw_score']=raw_ranks_df[0]+raw_ranks_df[1]
results['filtered_score']=filtered_ranks_df[0]+filtered_ranks_df[1]

results= pd.merge(results, drugs,how='left',on='_id')
columns_to_keep=['name','filtered_score','raw_score']
results=results[columns_to_keep]
results=results.sort_values(by='filtered_score',ascending=False)
results.to_csv('ComplEx_resutls.csv',index=False)
