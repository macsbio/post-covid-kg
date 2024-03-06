# coding: utf-8

"""Python file for querying the Disgenet database locally when the server is down"""

# Set path
import os 
from pathlib import Path
path = Path(__file__).resolve().parent
os.chdir(path)


from typing import Tuple
import pandas as pd
from pyBiodatafuse.utils import collapse_data_sources, get_identifier_of_interest


def get_disgenet_diseases(bridgedb_df: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    """Get location of gene in human body.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :returns: a DataFrame containing the OpenTargets output and dictionary of the query metadata.
    """
    data_df = get_identifier_of_interest(bridgedb_df, "Ensembl")
    data_df.drop(columns=['identifier.source','target.source'], inplace=True)
    
    data_conversion = get_identifier_of_interest(bridgedb_df, "NCBI Gene")
    data_conversion = data_conversion.rename(columns={'target': 'ncbi'})
    data_conversion= pd.merge(data_conversion,data_df,on='identifier',how='left') #add ensembl
    
    
    gene_ids_conversion = data_conversion[["ncbi",'target']]
    ###########################
    gene_ids_csv = pd.read_csv('disgenet/genes_disgenet.csv')
    diseases_ids_csv = pd.read_csv('disgenet/diseases_disgenet.csv')
    gene_diseases = pd.read_csv('disgenet/gene-disease_disgenet.csv')
    '''CURATED: GDAs from UniProt, PsyGeNET, Orphanet, the CGI, CTD (human data), ClinGen, and the Genomics England PanelApp.'''
    
    gene_diseases = gene_diseases[gene_diseases['source'].isin(['UNIPROT', 'PSYGENET', 'ORPHANET','CGI','CTD_human','CLINGEN','GENOMICS_ENGLAND'])]
    #str 
    gene_ids_csv['geneNID']=gene_ids_csv['geneNID'].astype(str)
    gene_ids_csv['geneId']=gene_ids_csv['geneId'].astype(str)
    gene_diseases['geneNID']=gene_diseases['geneNID'].astype(str)
    gene_diseases['diseaseNID']=gene_diseases['diseaseNID'].astype(str)
    diseases_ids_csv['diseaseNID']=diseases_ids_csv['diseaseNID'].astype(str)
    
    #Mapping
    df_all= pd.DataFrame()
    df_all= gene_ids_conversion
    df_all = df_all.rename(columns={'ncbi': 'geneId'})
    df_all = pd.merge(df_all,gene_ids_csv,on='geneId',how='left')
    df_all = pd.merge(df_all,gene_diseases,on='geneNID',how='left')
    df_all = pd.merge(df_all,diseases_ids_csv,on='diseaseNID',how='left')
    
    df_all=df_all[['diseaseId','geneId','diseaseName','score','geneName','target']]
    df_all = df_all.rename(columns={'diseaseName': 'disease_name'})
    df_all = df_all.rename(columns={'diseaseId': 'diseaseid'})
    
    df_all['geneId']=df_all['geneId'].astype(str)
    
    
    disgenet_df=df_all
    data_df = get_identifier_of_interest(bridgedb_df, "Ensembl")
    
    

    
    ##########################

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace="Ensembl",
        target_df=disgenet_df,
        common_cols=["target"],
        target_specific_cols=["disease_name", "diseaseid"],
        col_name="DisGeNET",
    )

    return merged_df

