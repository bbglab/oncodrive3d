import csv
import re
import numpy as np
import pandas as pd
import requests


## Parsers

def parse_maf_input(maf_input_path, keep_samples_id=False):
    """
    Parse in.maf file which is used as input for the HotMAPS method.
    """

    # Load
    maf = pd.read_csv(maf_input_path, sep="\t", dtype={'Chromosome': str})

    # Select only missense mutation and extract Gene name and mut
    maf = maf[maf['Variant_Classification'].str.contains('Missense_Mutation')]
    maf["Pos"] = maf.loc[:, "HGVSp_Short"].apply(lambda x: int(re.sub("\\D", "", (x[2:]))))
    maf["WT"] = maf["HGVSp_Short"].apply(lambda x: re.findall("\\D", x[2:])[0])
    maf["Mut"] = maf["HGVSp_Short"].apply(lambda x: re.findall("\\D", x[2:])[1])
    maf = maf[["Hugo_Symbol", "Pos", "WT", "Mut", "Tumor_Sample_Barcode"]]
    maf = maf.sort_values("Pos").rename(columns={"Hugo_Symbol" : "Gene"})

    if keep_samples_id == False:
        maf = maf.drop(columns=["Tumor_Sample_Barcode"])
    
    return maf.reset_index(drop=True)


def parse_cluster_output(out_cluster_path):
    """
    Parse out.gz file which is the output of HotMAPS.
    It will be used as ground truth for the new clustering method.
    """
    
    # Select only necessary info and rename to match the input file
    cluster = pd.read_csv(out_cluster_path, sep="\t")
    cluster = cluster.copy().rename(columns={"CRAVAT Res" : "Pos"})[["Pos", "Ref AA", "HUGO Symbol", 
                                                                     "Min p-value", "q-value"]]
    cluster = cluster.rename(columns={"HUGO Symbol" : "Gene", "Ref AA" : "WT"})
    
    return cluster


## Handle proteins fragments         ## IN THEORY DEPRECATED (TO CHECK)

def rounded_up_division(num, den):
    """
    Simply round up the result of the division.
    """
    
    return -(-num // den)


def get_pos_fragments(mut_gene_df):
    """
    Get the corresponding fragment of each position of the protein.
    """
    
    max_f = rounded_up_division(max(mut_gene_df["Pos"]), 1400)
    bins = [n * 1400 for n in range(0, max_f+1)]
    group_names = list(range(1, max_f+1))
    
    return pd.cut(mut_gene_df["Pos"], bins, labels = group_names)


## GENE-Uniprot mapping

def uniprot_to_hugo(uniprot_ids=None, hugo_as_keys=False, get_ensembl_id=False):
    """
    Returns a dictionary with UniProt IDs as keys and corresponding HUGO symbols as values.
    """
    # Download the latest HGNC gene symbols file
    url = "https://www.genenames.org/cgi-bin/download/custom?col=gd_app_sym&col=gd_pub_ensembl_id&col=md_prot_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"
    response = requests.get(url)
    csv_data = response.content.decode('utf-8')
    
    # Parse the CSV data into a dictionary of UniProt IDs to HUGO symbols
    reader = csv.reader(csv_data.splitlines(), delimiter='\t')
    next(reader)  # Skip the header row
    dict_output = {}
    for row in reader:
        hugo_symbol, ensembl_id, uniprot_id = row
        if uniprot_id and hugo_symbol:
            
            # Get Uniprot ID, HUGO symbol, Ensembl ID
            if get_ensembl_id:
                if not ensembl_id:
                    ensembl_id = np.nan
                
                if hugo_as_keys:
                    dict_output[hugo_symbol] = uniprot_id, ensembl_id
                else:
                    dict_output[uniprot_id] = hugo_symbol, ensembl_id
                    
            # Get Uniprot ID, HUGO symbol
            else:
                if hugo_as_keys:
                    dict_output[hugo_symbol] = uniprot_id
                else:
                    dict_output[uniprot_id] = hugo_symbol
    
    # Filter the dictionary to only include UniProt IDs that were input
    if uniprot_ids is not None and not hugo_as_keys:
        dict_output = {k: v for k, v in dict_output.items() if k in uniprot_ids}
    
    return dict_output


## Model confidence and PAE

def weighted_avg_plddt_vol(target_pos, mut_plddt_df, cmap):
    """
    Get the weighted average pLDDT score across residues in 
    the volume, based on number of mutations hitting each residue.  
    """
    
    return mut_plddt_df[[pos in np.where(cmap[target_pos-1])[0]+1 for pos in mut_plddt_df.Pos]].Confidence.mean()


def weighted_avg_pae_vol(target_pos, mut_plddt_df, cmap, pae):
    """
    Get the weighted average PAE across residues in the volume, 
    based on number of mutations hitting each residue.  
    """
    
    contacts = mut_plddt_df[[pos in np.where(cmap[target_pos-1])[0]+1 for pos in mut_plddt_df.Pos]].Pos.values
    pae_vol = pae[np.repeat(target_pos-1, len(contacts)), contacts-1].mean()
    
    return pae_vol


## Other utils

def get_samples_info(mut_gene_df, cmap):
    """
    Get total samples and ratio of unique samples having
    mutations in the volume of each mutated residues.
    """
    
    # Get total samples and # mutated samples of each mutated res
    tot_samples = len(mut_gene_df["Tumor_Sample_Barcode"].unique())
    pos_barcodes = mut_gene_df.groupby("Pos").apply(lambda x: x["Tumor_Sample_Barcode"].unique())
    pos_barcodes = pos_barcodes.reset_index().rename(columns={0 : "Barcode"})

    # Get the ratio of unique samples with mut in the vol of each mutated res
    uniq_pos_barcodes = [len(pos_barcodes[[pos in np.where(cmap[i-1])[0]+1 for 
                                           pos in pos_barcodes.Pos]].Barcode.explode().unique()) for i in pos_barcodes.Pos]
    pos_barcodes["Tot_samples"] = tot_samples        
    pos_barcodes["Samples_in_vol"] = uniq_pos_barcodes
    #pos_barcodes["Ratio_samples_in_vol"] = np.array(uniq_pos_barcodes) / tot_samples       
    
    return pos_barcodes


def get_unique_pos_in_contact(lst_pos, cmap):
    """
    Given a list of position and a contact map, return a numpy 
    array of unique positions in contact with the given ones.
    """
    
    return np.unique(np.concatenate([np.where(cmap[pos-1])[0]+1 for pos in lst_pos]))


def add_samples_info(mut_gene_df, result_pos_df, samples_info, cmap, pae=None):
    """
    Add information about the ratio of unique samples in the volume of 
    each mutated residues and in each detected community (meta-cluster) 
    to the residues-level output of the tool.
    """

    # Add samples in vol
    result_pos_df = result_pos_df.merge(samples_info.drop(columns=["Barcode"]), on="Pos", how="outer")
    
    # Get per-community ratio of mutated samples
    if result_pos_df["Cluster"].isna().all():
        result_pos_df["Samples_in_cl_vol"] = np.nan
        result_pos_df["Mut_in_cl_vol"] = np.nan
        result_pos_df["Res_in_cl"] = np.nan
        result_pos_df["pLDDT_cl_vol"] = np.nan
        #result_pos_df["Ratio_samples_in_comm"] = np.nan
    else:       
        community_pos = result_pos_df.groupby("Cluster").apply(lambda x: x.Pos.values)
        community_mut = community_pos.apply(lambda x: sum([pos in get_unique_pos_in_contact(x, cmap) for 
                                                           pos in mut_gene_df.Pos]))
        community_samples = community_pos.apply(lambda x: 
                                        len(mut_gene_df[[pos in get_unique_pos_in_contact(x, cmap) for 
                                                         pos in mut_gene_df.Pos]].Tumor_Sample_Barcode.unique()))
        community_plddt = community_pos.apply(lambda x: mut_gene_df.Confidence[[pos in get_unique_pos_in_contact(x, cmap) 
                                                                             for pos in mut_gene_df.Pos]].mean())

        #community_samples_ratio = community_samples / samples_info["Tot_samples"].unique()[0]
        #community_mut_ratio = community_mut / len(mut_gene_df)
        community_pos_count = community_pos.apply(lambda x: len(x))
        community_samples = pd.DataFrame({"Samples_in_cl_vol" : community_samples, 
                                          #"Ratio_samples_in_comm" : community_samples_ratio,
                                          "Mut_in_cl_vol" : community_mut,
                                          #"Ratio_mut_in_comm" : community_mut_ratio,
                                          "Res_in_cl" : community_pos_count,
                                          "pLDDT_cl_vol" : community_plddt})
        
        # Add to residues-level result
        result_pos_df = result_pos_df.merge(community_samples, on="Cluster", how="outer")
        
    # AF PAE
    if pae is not None:
        result_pos_df["PAE_vol"] = result_pos_df.apply(lambda x: weighted_avg_pae_vol(x["Pos"], mut_gene_df, cmap, pae), axis=1)  
    else:
        result_pos_df["PAE_vol"] = np.nan
        
    # AF confidence
    result_pos_df["pLDDT_res"] = result_pos_df.apply(lambda x: mut_gene_df.Confidence[mut_gene_df["Pos"] == x.Pos].values[0], axis=1)
    result_pos_df["pLDDT_vol"] = result_pos_df.apply(lambda x: weighted_avg_plddt_vol(x["Pos"], mut_gene_df, cmap), axis=1)
    result_pos_df["pLDDT_cl_vol"] = result_pos_df.pop("pLDDT_cl_vol")

    # Sort positions
    result_pos_df = result_pos_df.sort_values("Rank").reset_index(drop=True)
    
    return result_pos_df


def add_nan_clust_cols(result_gene):
    """
    Add columns showing clustering results with only nan for 
    genes that are not tested (not enough mutations, etc).
    """

    result_gene = result_gene.copy()
    
    columns = ["pval", "qval", "C_gene", "C_pos", 'C_label', 'Ratio_obs_sim_top_vol', 
               "Clust_res", 'Clust_mut', 'Mut_in_top_vol', "Mut_in_top_cl_vol",
               'Tot_samples', 'Samples_in_top_vol', 'Samples_in_top_cl_vol', 
               "PAE_top_vol", "pLDDT_top_vol", "pLDDT_top_cl_vol", 'F']
    
    for col in columns:
        result_gene[col] = np.nan

    return result_gene


def sort_cols(result_gene):
    """
    Simply change the order of columns of the genes-level result.
    """

    cols = ['Gene', 'Uniprot_ID', 
            'pval', 'qval', 'C_gene', 'C_pos', 'C_label', 'Ratio_obs_sim_top_vol', "Clust_res",
            'Mut_in_gene', 'Clust_mut', 'Mut_in_top_vol', "Mut_in_top_cl_vol",
            'Tot_samples', 'Samples_in_top_vol', 'Samples_in_top_cl_vol', 
            "PAE_top_vol", "pLDDT_top_vol", "pLDDT_top_cl_vol",
            'F', 'Status', 'Cancer', 'Cohort']

    return result_gene[[col for col in cols if col in result_gene.columns]]