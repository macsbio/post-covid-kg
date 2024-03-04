#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 12:50:01 2024

@author: alejandroadriaquelozano
"""


import requests
import pandas as pd
from typing import Tuple
import math

def get_drug_disease_interactions(drug_df: pd.DataFrame,disgenet_result:pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    base_url = "https://api.platform.opentargets.org/api/v4/graphql"
    """Get information about drugs associated with diseases of interest.

    :param drug:df: get_gene_disease_associations output for creating the list of diseases ids to query
    :returns: a DataFrame containing the 
    """

    # Iterate through the dictionary and remove entries with NaN values for drugs and diseases
    #Drugs
    drugs_all = dict(drug_df['ChEMBL_Drugs'])
    drugs_list=[]
    for key, value in drugs_all.items():
        for lists in value:
             drugs_list.append(lists['chembl_id'])

    drugs_ids = [value for value in set(drugs_list) if not (math.isnan(value) if isinstance(value, float) else False)]
  
    #Diseases
    diseases_all = dict(disgenet_result['DisGeNET'])
    diseases_list=[]
    for key, value in diseases_all.items():
        for lists in value:
             diseases_list.append(lists['diseaseid'])

    diseases_ids = [value for value in set(diseases_list) if not (math.isnan(value) if isinstance(value, float) else False)]

   
    diseases_drugs_list= []
    for x in drugs_ids:
        query_string = """
                        query KnownDrugsQuery(
              $cursor: String
              $freeTextQuery: String
              $size: Int = 10
            ) {
              drug(chemblId: $chemblId) {
                id
                name
                knownDrugs(cursor: $cursor, freeTextQuery: $freeTextQuery, size: $size) {
                  count
                  rows {
                    disease {
                      id
                      name
                      dbXRefs
                    }
                    target {
                      id
                      approvedName
                      approvedSymbol
                    }
                  }
                }
              }
            } 
            
            """
        query_string = query_string.replace("$chemblId",  '"' + x + '"')
        r = requests.post(base_url, json={"query": query_string}).json()
        diseases_drugs_list.append(r['data'])
    
    
    # Extracting and simplifying the structure
    col_names = ['identifier','drug_name','drug_diseases']
    drug_disease_df = pd.DataFrame(columns=col_names)
    
        # Extracting and simplifying the structure
    col_names = ['identifier', 'drug_name', 'drug_diseases']
    drug_disease_df = pd.DataFrame(columns=col_names)
    
    for entry in diseases_drugs_list:
        entries = entry['drug']
        drugs = entries['knownDrugs']['rows']
        drugs_list = []
        dict_new = None  # Initialize dict_new outside the loop
        
        for entry in drugs:
            drugs_fixed = entry['disease']
            other_IDs = drugs_fixed['dbXRefs']
            
            umls_id = None  # Initialize umls_id before the loop
            
            for value in other_IDs:
                if 'UMLS' in value:
                    result_list = value.split(':')
                    umls_id = result_list[1]
                    break  # Once we find the UMLS ID, we can exit the loop
                    
            if umls_id is None:
                umls_id = math.nan
            
            del drugs_fixed['dbXRefs']
            drugs_fixed['umls'] = umls_id
            
            if drugs_fixed['umls'] in diseases_ids:  # only include diseases in my graph
                drugs_list.append(drugs_fixed)
        
        if drugs_list:
            dict_new = {'identifier': entries['id'], 'drug_name': entries['name'], 'drug_diseases': drugs_list}
        
        if dict_new is not None:
            dict_new_df = pd.DataFrame([dict_new])
            drug_disease_df = pd.concat([drug_disease_df, dict_new_df], ignore_index=True)
    
    return drug_disease_df
