#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 16:30:05 2024

@author: alejandroadriaquelozano
"""
# Set path
import os 
from pathlib import Path
path = Path(__file__).resolve().parent
os.chdir(path)

# Import modules
import pickle
import pandas as pd
from pyBiodatafuse import id_mapper
from pyBiodatafuse.annotators import opentargets, stringdb, wikipathways
from pyBiodatafuse.utils import combine_sources
import generator,neo4j_exporter, new_disgenet_annotator,drug_disease_annotator,minerva

# Differentially expressed genes in post-COVID-19
genes_of_interest = """LOC729609
LOC105374060
DMP1
PNLIP
OR4N3P
SLC6A14
LOC101927239
DEFB105A
DEFB105B
GSTTP1
NEUROD1
RND1
VN1R10P
LOC440446
LOC152225
LOC101929341
PGLYRP3
LINC01533
LINC01090
SPEM1
C16orf82
MIR4432HG
LINC01169
FAM71A
RNASE10
KLF17
C9
ARC
MYL10
GCM1
AIPL1
HSPA6
LOC101929124
C7orf65
SLC2A14
PNLIPRP2
NPAS4
LOC101060498
PROP1
ELAVL3
LOC105747689
TNF
ADAMTS4
PCDH10
LOC101927274
NR4A2
LOC102724612
CEACAM22P
SNAI1
SLC2A3
DLX3
ID2
LOC151475
ATF3
NKAIN4
ASAP1-IT2
NOXRED1
DNM1P41
SLC7A11
C10orf82
ULBP2
TPTE2P6
NR4A3
LOC399715
CNTN3
GEM
HSPA7
NCMAP
PNP
PLK2
ATP2C2
TNFRSF10D
ULBP3
HSPA5
EFHB
HSD17B13
WNK3
LINC01535
ELL2
RND3
DUSP5
NRXN3
IPCEF1
ZNF492
SDR16C5
CENPL
SOX11
MAFF
PRG4
PCDH17
CDKN1A
PELI1
TMEM169
TMEM236
EFNA5
GCH1
ANGPTL4
MAP1LC3C
CHL1
MPZ
SERPINE1
SLC2A1
LRRC16A
FRZB
GLIS3
TIAM1
SRGAP1
SH2D4A
MYEF2
NT5E
VGLL3
PRTG
DPP4
KLF11
TAF13
STRADB
POMP
LAMTOR5
CCDC69
ZNF32
IQSEC2
APIP
GDF9
SCUBE2
C20orf24
ZSWIM7
TIMM8B
LOC102724532
PRR16
AHRR
LEFTY2
IRX3
VMO1
PVALB
MT1DP
CALML5
LOC101929116
LOC101929694
LINC01205
LINC01241
TMPRSS11A
LOC101928942
LOC100507461
LINC01565
LOC101928358
SCGB1D4
TTR
LINC01284
SSX8
TMEM225
NCRNA00250
OR13D1
LINC01192
CALCB
LINC00411
LINC01227
MIR5689HG
LINC00615
GHSR
LOC105375556
CT45A5
LOC646029
ZFP42
CT45A9
FLJ46066
CGA
LOC285692
LOC105369509
CLEC1B
HIST1H4A
DSCAM-IT1
CT45A2
CT45A8
LINC00928
BDKRB1
LOC105370586
TRIM51
LOC101927480
LINC01568
CASC17
LOC101929631
LINC01233
LOC101927948
OR13C5
SSX2
SSX2B
CACNA1C-IT3
LOC100500773
SPATA3
LOC101927374
FBXO47
LINC01493
LOC105369431
LOC105376468
OR5W2
REG4
CD5L
LINC01514
LOC105376331
LOC102723557
PISRT1
HIGD2B
PAGE1
MMP26
LOC101928602
LOC102723895
ACTR3BP2
LOC101927363
HNRNPKP3
LOC101927188
DISC1-IT1
LOC102467222
FAM9B
GLOD5
C2orf48
LOC100288254
FRG2
GACAT3
FOXCUT
LOC101927357
LOC101929260
OR13C2
LOC101929754
LOC146513
OR2AT4
PBOV1
TFDP3
LOC101929420
HRAT17
OR6W1P
SSX9
SSX3
HMGA1P7
LINC00374
LINC01288
LINC00836
LINC01320
TRIM64
SDR16C6P
LOC729966
LOC105375014
LINC01441
SCNN1G
C7orf69
OPN1LW
KRTAP5-4
ANKUB1
TMEM213
TFAP2D
DANT2
LOC101927419
TXNDC2
OR11A1
LINC01317
LOC101805491
LOC286083
LOC101929563
LINC01216
LINC01163
LOC101927166
PHOX2B
LOC102467081
CT45A6
SND1-IT1
SSX4B
SSX4
SULT1E1
NOL4
ZNF716
SUMO1P1
LOC440896
G6PC
MIR31HG
LOC101929259
HTR3C
LOC730100
MAB21L3
IL6
MIP
TRIM64B
CNGB1
LINC01531
FOXL2NB
CXCL8
SLC15A1
GABRB1
LINC00862
ZPBP2
LOC101928992
DPPA4
POU2F3
NUTM1
LOC105372440
SELE
GPR143
FSTL5
AXDND1
LINC01619
SAMD7
LOC100131257
ABCC13
C17orf78
CRX
C12orf42
FOXG1
HTR3A
LOC644189
PNPLA1
LINC00880
TOP1P2
CAGE1
LINC00670
LOC101928231
FAM138C
RTP1
LOC101928617
SPAG11B
LOC101927691
SLC35G3
BCO1
SLC35G4
LINC00636
EPGN
PTGS2
PGC
LOC102724467
LOC101928103
TRPC5OS
LOC338694
LINC01036
DLX6
LINC00426
CXorf65
HP09025
LOC389273
DPCR1
C5orf60
PCSK1
LOC494141
GADD45B
C1orf87
ANKS4B
JAKMIP2
LINC00266-3
DRAIC
TCAM1P
MIR202HG
SPRR2F
FAM138B
LINC00907
CCL19
ASCL1
NUP210L
LINC01170
LINC00264
ANKRD7
LOC102724601
SH2D6
FAM138F
FAM138A
GYPE
DDX4
IL5RA
TNFRSF9
LINC00368
LGSN
NEK5
LOC105374177
GLB1L3
LOC105379511
MT1A
FAM138E
TEKT3
SV2C
NR2E3
PLA2G10
LOC101927770
ENO4
SBK2
A2ML1
LOC101927257
SPRY4-IT1
DNAH8
AK7
ASXL3
TEX38
DNM1P35
CCL26
PPP3R2
CTSLP2
ACBD7
SOX2-OT
STC1
LOC284865
FDPSP2
MARVELD2
CDKL2
DCX
SHISA9
C4orf26
DNAH5
CD3G
TTC23L
PDE6A
APOBEC3H
LINC00311
CXCL2
LINC00632
SALL4
LOC105372582
FAM106CP
RASD1
CACNA1F
ELAVL2
KIAA0087
GIPR
CIDEA
BCL11B
TNFRSF11B
CA13
ANKRD20A9P
FAM106B
SEMA3E
GPRC5A
LOC285819
LOC730101
IL1RL1
RGS2
RYBP
C3orf52
HOOK1
PCDH9
CDH19
PGA4
STARD4
CYP2B7P
TFPI2
PDK4
PGA5
KCNAB3
LINC00641
LOC102724571
SEZ6L
TNFSF9
ZNF483
M1AP
FAAP24
KLHL15
CHD1
AP1S3
CDS1
CRTAC1
GYG2
GRHL1
FSIP1
SYT1
PLCXD3
LOC101928371
PEG10
MPZL3
ZNF331
KCNQ1OT1
LOC388436
LOC79999
FAM106A
RPS6KA6
BCL2L15
TBX5
EMP1
PPP2R2B
TACR1
SLC7A10
ELOVL6
ATP1B3
SEMA4A
CEP152
LINC01296
NRXN1
ADGRG2
CLDN1
ZSWIM6
WNT3
CCDC170
THBS1
SLC35F2
ZC3H12B
PLIN1
LOC401052
CATSPERG
IFRD1
GAS2L3
APOBEC3D
POU2F2
ERRFI1
ARSJ
FOXC1
PRDM1
RASGRP1
KIAA1683
PRELP
TIPARP
ZC3H12A
SGIP1
PDE8B
GFPT2
CABP4
RAD51B
MICB
EIF4A3
FAM72C
C7
QPCT
MAP3K8
TUFT1
DUXAP10
SHROOM3
ZC3HAV1
S1PR2
FAM122C
HRH1
UGCG
SOX9
LYVE1
BCL2L11
EIF2AK3
C11orf63
SERPINB8
LEPR
CACNB2
CACNA2D4
NR2F1
CLCF1
PSD3
ADNP2
DYNC2H1
OR2A20P
SYT17
VASH2
TMEM2
OR2A9P
USP32P2
EDIL3
LOX
MXD1
NHSL1
DLC1
CYBB
ETV5
CEP126
PTPRF
COCH
SCRN1
PPM1D
LILRB4
MFSD4A
CCDC144B
PXDNL
AHR
TRIM14
FRMD4B
CD84
TIAM2
ADAMTS5
XYLT1
MYOF
SLC7A1
SMG1P3
UGDH
PMP22
AMPH
NPIPB5
NT5DC3
UBE2D2
PIGX
TTC1
SRP14
GKAP1
FIBP
MED11
VTI1B
ATPAF1
DNAJC19
MRPL24
TRIM16L
POLR2F
GCSH
TMEM147
LSM10
MRPL40
C11orf74
SERF2-C15ORF63
NDUFAF2
UBE3D
MALSU1
COA4
ELP6
MTX2
CMC4
MON1A
CABP7
MID1IP1
COA6
KIF22
TSEN15
NDFIP2
HYPK
ZCRB1
PARK7
COX16
GTF3C6
MINOS1
MRPS15
STOML2
KCNS3
CACNA2D3
CTNNBIP1
C7orf55
COPS5
CHCHD5
YBX3P1
SPAG7
NDUFS3
TPI1
PET100
ST3GAL2
MRPL21
TP53TG1
CDKN2AIPNL
OIP5
RPS20
ATP5E
CBWD2
CDK5
TOMM5
PRR34
HINT1
BAD
ATP5L
SFXN5
AAMDC
MRPL51
KIAA0930
VAMP5
SEPW1
NDUFA6
SLIRP
SHISA2
NUDT2
COX5B
SNRPN
SNURF
AURKA
CBWD1
NDUFB2
NAA38
CKM
GPD1
RPS29
DHRS4L1
MRPL33
LOC100507291
ATP23
UQCRQ
NDUFC2
BOLA3
TCEB2
COX7A1
DHRS4
COX6C
FHL2
SLN
NDUFA1
RPL21P28
RPL21
NDUFC2-KCTD14
ATP5I
UQCC2
LOC101929231
DBNDD1
NDUFB9
LAMB3
CSF3R
USMG5
DHRS4L2
SERPINA1
C1orf53
GLT1D1
GREM2
UQCRBP1
FAM24B
S100A8
CDH22
LEFTY1"""

gene_list = genes_of_interest.split("\n")
len(gene_list)

data_input = pd.DataFrame(gene_list, columns=["identifier"])
data_input.head()

bridgdb_df, bridgdb_metadata = id_mapper.bridgedb_xref(
    identifiers=data_input,
    input_species="Human",
    input_datasource="HGNC",
    output_datasource="All",
)
bridgdb_df.head()


disgenet_result = new_disgenet_annotator.get_disgenet_diseases(bridgdb_df)

loc_df, opentargets_loc_metadata = opentargets.get_gene_location(bridgedb_df=bridgdb_df)
loc_df.head()

go_process_df, opentargets_go_metadata = opentargets.get_gene_go_process(bridgedb_df=bridgdb_df)
go_process_df.head()

reactome_process_df, opentargets_process_metadata = opentargets.get_gene_reactome_pathways(
    bridgedb_df=bridgdb_df
)
reactome_process_df.head()

drug_df, opentargets_drug_metadata = opentargets.get_gene_drug_interactions(bridgedb_df=bridgdb_df)
drug_df.head()

disease_df, opentargets_disease_metadata = opentargets.get_gene_disease_associations(
    bridgedb_df=bridgdb_df
)
disease_df.head()

wp_df, wp_metadata = wikipathways.get_gene_wikipathway(bridgedb_df=bridgdb_df)
wp_df.head()

ppi_df, ppi_metadata = stringdb.get_ppi(bridgedb_df=bridgdb_df)
ppi_df.head()

#MINERVA
project_list_df = minerva.list_projects()
map_components = minerva.get_minerva_components(project_list_df, map_name='COVID19 Disease Map', get_reactions=False )
minerva_df = minerva.get_gene_minerva_pathways(bridgdb_df,map_components)

# DRUG DISEASE RELATIONSHIPS
drug_disease = drug_disease_annotator.get_drug_disease_interactions(drug_df,disgenet_result)


combined_df = combine_sources(
    [
        disgenet_result,
        loc_df,
        go_process_df,
        reactome_process_df,
        drug_df,
        disease_df,
        wp_df,
        ppi_df,
        minerva_df
    ]
)

combined_df.to_csv('graph.csv')

with open("combined_df.pkl", "wb") as out:
    pickle.dump(combined_df, out)
    
combined_df = generator.load_dataframe_from_pickle("combined_df.pkl")

pygraph = generator.generate_networkx_graph(combined_df,drug_disease)

neo4j_exporter.save_graph_to_neo4j_graphml(pygraph, 'graph.graphml')
