import json
import os
import time

import requests
import numpy as np
import daiquiri
import pandas as pd
from tqdm import tqdm
from scripts import __logger_name__

logger = daiquiri.getLogger(__logger_name__ + ".plotting.uniprot_feat")


def get_evidence(feat):
    """
    Get source of evidence ID and reference.
    """

    if "evidences" in feat.keys():
        evidence = list(set([f'{e["source"]["name"] if "source" in e.keys() else np.nan}' for e in feat["evidences"]]))
    else:
        evidence = np.nan
    return evidence


def get_domain_id(feat):
    """
    Get domain ID.
    """

    if "evidences" in feat.keys() and "source" in feat["evidences"][0]:
        domain_id = feat["evidences"][0]["source"]["id"]
    else:
        domain_id = np.nan
    return domain_id


def get_description(feat):
    """
    Get feature description.
    """

    if "description" in feat.keys():
        description = feat["description"]
        if len(description) == 0:
            description = np.nan
    else:
        description = np.nan
    return description


def _uniprot_request_feat(lst_uniprot_ids):
    """
    Use Features from EMBL-EBI Proteins API to get 
    a json including protein features.
    
    https://www.ebi.ac.uk/proteins/api/doc/#featuresApi
    https://doi.org/10.1093/nar/gkx237
    """
    
    prot_request = [f"{prot}" if i == 0 else f"%2C{prot}" for i, prot in enumerate(lst_uniprot_ids)]
    requestURL = f"https://www.ebi.ac.uk/proteins/api/features?offset=0&size=100&accession={''.join(prot_request)}"
    
    status = "INIT"
    while status != "FINISHED":
        if status != "INIT":
            time.sleep(10)
        try:
            r = requests.get(requestURL, headers={ "Accept" : "application/json"}, timeout=160)
            if r.ok:
                status = "FINISHED"
            else:
                logger.debug(f"Error occurred after successfully sending request. Status: {r.raise_for_status()}")             
                status = "ERROR"
        except requests.exceptions.RequestException as e:                          
            status = "ERROR"                                                     
            logger.debug(f"Request failed: {e}")                               
    
    for dictio in json.loads(r.text):

        yield dictio
        

def get_batch_prot_feat(batch_ids):
    """
    Parse the json obtained from the Features 
    service extracting protein features.
    
    https://www.ebi.ac.uk/proteins/api/doc/#featuresApi
    https://doi.org/10.1093/nar/gkx237
    """
    
    lst_uni_id = []
    lst_type = []
    lst_begin = []
    lst_end = []
    lst_description = []
    lst_evidence = []
    lst_domain_id = []

    types = ["DOMAIN", "DNA_BIND", 
             "ACT_SITE", "BINDING", "SITE",
             "MOD_RES", "CARBOHY", "LIPID", "CARBOHYD", "CROSSLNK",
             "MOTIF", "ZN_FING", 'TRANSMEM', 'INTRAMEM', 'SIGNAL']
    
    for protein in _uniprot_request_feat(batch_ids):
        uni_id = protein["accession"]

        for feat in protein["features"]:
            if feat["type"] in types:
                lst_uni_id.append(uni_id)
                lst_type.append(feat["type"])
                lst_begin.append(feat["begin"])
                lst_end.append(feat["end"])
                lst_description.append(get_description(feat))
                lst_evidence.append(get_evidence(feat))
                if feat["type"] == "DOMAIN":
                    lst_domain_id.append(get_domain_id(feat))
                else:
                    lst_domain_id.append(np.nan)
                    
    return pd.DataFrame({"Uniprot_ID" : lst_uni_id, 
                         "Type" : lst_type, 
                         "Begin" : lst_begin, 
                         "End" : lst_end, 
                         "Description" : lst_description, 
                         "Evidence" : lst_evidence, 
                         "Domain_ID" : lst_domain_id})


def get_prot_feat(ids, batch_size=100):
    """
    Use the Features service from Proteins API of EMBL-EBI to get 
    protein features of all provided Uniprot IDs.
    
    https://www.ebi.ac.uk/proteins/api/doc/#featuresApi
    https://doi.org/10.1093/nar/gkx237
    """
    
    lst_df = []
    batches_ids = [ids[i:i+batch_size] for i in range(0, len(ids), batch_size)]

    for batch_ids in tqdm(batches_ids, total=len(batches_ids), desc="Extracting protein features from Uniprot"):
        
        batch_df = get_batch_prot_feat(batch_ids)
        lst_df.append(batch_df)
    
    return pd.concat(lst_df).reset_index(drop=True)


def parse_prot_feat(feat_df):
    """
    Parse dataframe including Features obtained by Proteins API of EMBL-EBI.
    Merge similar entries to simplify visualization of the result. 

    https://www.ebi.ac.uk/proteins/api/doc/#featuresApi
    https://doi.org/10.1093/nar/gkx237
    """

    feat_df = feat_df.copy()
        
    # Add PTM description for PTM
    feat_df["Full_description"] = np.nan
    feat_df["Full_description"] = feat_df["Description"]
    feat_df["Evidence"] = feat_df.pop("Evidence")

    # PTMs
    phosp_ix = feat_df['Description'].str.contains('Phosp', case=False).fillna(False)
    acetyl_ix = feat_df['Description'].str.contains('Acetyl', case=False).fillna(False)
    methyl_ix = feat_df['Description'].str.contains('Methyl', case=False).fillna(False)
    feat_df.loc[(feat_df["Type"] == "MOD_RES") & phosp_ix, "Description"] = "Phosphorilation"
    feat_df.loc[(feat_df["Type"] == "MOD_RES") & acetyl_ix, "Description"] = "Acetylation"
    feat_df.loc[(feat_df["Type"] == "MOD_RES") & methyl_ix, "Description"] = "Methylation"
    feat_df.loc[(feat_df["Type"] == "MOD_RES") & ~phosp_ix & ~acetyl_ix & ~methyl_ix, "Description"] = "Others"

    # Other PTMs
    feat_df.loc[feat_df["Type"] == "MOD_RES", "Description"] = feat_df[feat_df["Type"] == "MOD_RES"].apply(lambda x: x["Description"].split(";")[0], axis=1)
    feat_df.loc[feat_df["Type"] == "MOD_RES", "Type"] = "PTM"
    feat_df.loc[feat_df["Type"] == "CARBOHYD", "Description"] = "Glycosylation"
    feat_df.loc[feat_df["Type"] == "CARBOHYD", "Type"] = "PTM"
    feat_df.loc[feat_df["Type"] == "LIPID", "Description"] = "Lipidation"
    feat_df.loc[feat_df["Type"] == "LIPID", "Type"] = "PTM"

    # Cross-links PTMs
    ubiqui_ix = feat_df['Description'].str.contains('Ubiquitin', case=False).fillna(False)
    sumo_ix = feat_df['Description'].str.contains('Sumo', case=False).fillna(False)
    feat_df.loc[(feat_df["Type"] == "CROSSLNK") & ubiqui_ix, "Description"] = "CL-Ubiquitination"
    feat_df.loc[(feat_df["Type"] == "CROSSLNK") & sumo_ix, "Description"] = "CL-SUMOylation"
    feat_df.loc[(feat_df["Type"] == "CROSSLNK") & ~ubiqui_ix & ~sumo_ix, "Description"] = "CL-Others"
    feat_df.loc[feat_df["Type"] == "CROSSLNK", "Type"] = "PTM"

    # Membrane
    feat_df.loc[feat_df["Type"] == "INTRAMEM", "Description"] = "Intra"
    feat_df.loc[feat_df["Type"] == "TRANSMEM", "Description"] = "Trans"
    feat_df.loc[(feat_df["Type"] == "INTRAMEM") | (feat_df["Type"] == "TRANSMEM"), "Type"] = "MEMBRANE"

    # Sites
    cleavage_ix = feat_df['Description'].str.contains('Cleavage', case=False).fillna(False)
    interaction_ix = feat_df['Description'].str.contains('Interaction', case=False).fillna(False)
    breakpoint_ix = feat_df['Description'].str.contains('Breakpoint', case=False).fillna(False)
    ubiquit_ix = feat_df['Description'].str.contains('Ubiquit', case=False).fillna(False)
    fusion_ix = feat_df['Description'].str.contains('Fusion point', case=False).fillna(False)

    feat_df.loc[(feat_df["Type"] == "SITE") & cleavage_ix, "Description"] = "Cleavage"
    feat_df.loc[(feat_df["Type"] == "SITE") & interaction_ix, "Description"] = "Interaction"
    feat_df.loc[(feat_df["Type"] == "SITE") & breakpoint_ix, "Description"] = "Breakpoint"
    feat_df.loc[(feat_df["Type"] == "SITE") & ubiquit_ix, "Description"] = "Ubiquitin"
    feat_df.loc[(feat_df["Type"] == "SITE") & fusion_ix, "Description"] = "Fusion point"
    feat_df.loc[(feat_df["Type"] == "SITE") & ~cleavage_ix & ~interaction_ix 
                & ~breakpoint_ix & ~ubiquit_ix & ~fusion_ix & ~cleavage_ix, "Description"] = "Others"
    
    feat_df.loc[feat_df["Type"] == "ACT_SITE", "Description"] = "Active"
    feat_df.loc[feat_df["Type"] == "BINDING", "Description"] = "Binding"
    feat_df.loc[(feat_df["Type"] == "ACT_SITE") | (feat_df["Type"] == "BINDING") | (feat_df["Type"] == "SITE"), "Type"] = "SITE"

    # Motifs    
    sumo_ix = feat_df['Description'].str.contains('Sumo', case=False).fillna(False)
    feat_df.loc[(feat_df["Type"] == "MOTIF") & sumo_ix, "Description"] = "SUMO-related"        
    feat_df.loc[(feat_df["Type"] == "MOTIF") & ~sumo_ix, "Description"] = "Others"        
    feat_df.loc[feat_df["Type"] == "ZN_FING", "Description"] = "Zinc finger"
    feat_df.loc[(feat_df["Type"] == "MOTIF") | (feat_df["Type"] == "ZN_FING"), "Type"] = "MOTIF"
        
    # Regions
    feat_df.loc[feat_df["Type"] == "SIGNAL", "Description"] = "Signal peptide"
    feat_df.loc[feat_df["Type"] == "DNA_BIND", "Description"] = "DNA binding"
    feat_df.loc[(feat_df["Type"] == "SIGNAL") | (feat_df["Type"] == "DNA_BIND"), "Type"] = "REGION"
    
    # Domain
    feat_df.loc[feat_df["Type"] == "DOMAIN", "Description"] = feat_df[feat_df["Type"] == "DOMAIN"].apply(
        lambda x: x["Description"].split(";")[0] if pd.notna(x["Description"]) else np.nan, axis=1)
    feat_df.loc[feat_df["Type"] == "DOMAIN", "Description"] = feat_df.loc[feat_df["Type"] == "DOMAIN", 
                                                                          "Description"].str.replace(r' \d+', '')
    feat_df["Domain_ID"] = feat_df.pop("Domain_ID")
    
    return feat_df


def add_feat_metadata(feat_df, seq_df):
    
    # Add metadata to Uniprot Feat
    feat_df["Evidence"] = feat_df["Evidence"].astype(str)
    feat_df = seq_df[["Gene", "Uniprot_ID", "Ens_Transcr_ID", "Ens_Gene_ID"]].merge(
        feat_df, how="left", on=["Uniprot_ID"]).drop_duplicates()
    feat_df = feat_df.dropna(how="any", subset=["Begin", "End"]).reset_index(drop=True)
    
    # Parse weird end positions
    feat_df = feat_df.copy() 
    feat_df = feat_df[feat_df["End"] != "~"]
    feat_df["End"] = feat_df["End"].astype(str).str.replace("~", "")
    feat_df["End"] = feat_df["End"].astype(str).str.replace(">", "")
    feat_df["Begin"] = feat_df["Begin"].astype(str).str.replace("<", "")
    feat_df["Begin"] = feat_df["Begin"].astype(str).str.replace("~", "")
    feat_df["Begin"] = pd.to_numeric(feat_df["Begin"], errors='coerce')
    feat_df["End"] = pd.to_numeric(feat_df["End"], errors='coerce')
    feat_df[["Begin", "End"]] = feat_df[["Begin", "End"]].astype(int)
    
    return feat_df


def get_uniprot_feat(seq_df, pfam_df, output_tsv):
    """
    Extract and parse dataframe including Features obtained by Proteins API of EMBL-EBI.
    Merge similar entries to simplify visualization of the result. 
    Add Pfam domain, HUGO symbol, Ensembl Gene and Transcript info.
    
    https://www.ebi.ac.uk/proteins/api/doc/#featuresApi
    https://doi.org/10.1093/nar/gkx237
    """
    
    feat_df = get_prot_feat(seq_df.Uniprot_ID)  
    feat_df = parse_prot_feat(feat_df)                     
    feat_df = add_feat_metadata(feat_df, seq_df)
    feat_df = pd.concat((feat_df, pfam_df)).sort_values(["Gene", "Uniprot_ID", "Begin"]).reset_index(drop=True)
    feat_df.to_csv(output_tsv, sep="\t", index=False)
    logger.debug(f"Uniprot Features are saved to {output_tsv}")