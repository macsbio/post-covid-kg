#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 15:53:47 2024

@author: alejandroadriaquelozano
"""
# Set path
import os 
from pathlib import Path
path = Path(__file__).resolve().parent
os.chdir(path)

import xml.etree.ElementTree as ET
from tqdm import tqdm
import pandas as pd

dB_file = 'full database drugbank.xml'
organism = 'Humans'
saveFile = 'atc_hierarchy.csv'


# Parse the XML file
xtree = ET.parse(dB_file)
xroot = xtree.getroot()

drugs = list(xroot)

df= pd.DataFrame()
for i in tqdm(range(len(drugs))):
    drug = drugs[i]
    idDB = drug[0].text # Drug Bank ID
    for idx,feature in enumerate(drug):
        if 'external-identifiers' in str(feature): #other drug's IDs
            aux = [ext_id[0].text + ":" + ext_id[1].text \
                                        for ext_id in list(drug[idx])]
                
            drug_external_id = ';'.join(aux)
        if 'atc-codes' in str(feature): # drug name
            atc_more = list()
            
            for atc in list(drug[idx]):
                code = atc.attrib['code']
                sub_codes=[]
                for i in range(len(atc)):
                    sub_code = atc[i].attrib['code']
                    sub_codes.append(sub_code)
                codes_list_entry = [code,*sub_codes,'Drug']      
            
                combined_codes = ','.join(codes_list_entry)
                atc_more.append(combined_codes)
            atc_more= ';'.join(atc_more)
            if len(atc_more) == 0 :
                 atc_more='Drug'
                
    row= {'dbID':idDB,'atc':atc_more,"ext_id":drug_external_id,}
    row_df= pd.DataFrame([row])
    df = pd.concat([df, row_df], ignore_index=True)




row=0
df['id']=None
for x in df['ext_id']:
    x=x.split(';')
    for theid in x:
        if 'ChEMBL' in theid:
            final= theid.split(':')
        
            df['id'][row]=final[1]
    row=row+1
        
df=df.drop(columns=['ext_id'])
df = df.rename(columns={'atc': 'atcClassification'})
        
df.to_csv(saveFile)
        