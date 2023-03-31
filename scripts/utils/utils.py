import re
import pandas as pd
import numpy as np
import csv
import requests


## Parsers

def parse_maf_input(maf_input_path, keep_samples_id=False):
    """
    Parse in.maf file which is used as 
    input for the HotMAPS method.
    """

    # Load
    maf = pd.read_csv(maf_input_path, sep="\t", dtype={'Chromosome': str})

    # Select only missense mutation and extract Gene name and mut
    maf = maf.loc[maf.Variant_Classification == "Missense_Mutation"].copy()
    maf["Pos"] = maf.loc[:, "HGVSp_Short"].apply(lambda x: int(re.sub("\\D", "", (x[2:]))))
    maf["WT"] = maf["HGVSp_Short"].apply(lambda x: re.findall("\\D", x[2:])[0])
    maf["Mut"] = maf["HGVSp_Short"].apply(lambda x: re.findall("\\D", x[2:])[1])
    maf = maf[["Hugo_Symbol", "Pos", "WT", "Mut", "Tumor_Sample_Barcode"]]
    maf = maf.sort_values("Pos").rename(columns={"Hugo_Symbol" : "Gene"})

    if keep_samples_id == False:
        maf = maf.drop(columns=["Tumor_Sample_Barcode"])
    
    return maf


def parse_cluster_output(out_cluster_path):
    """
    Parse out.gz file which is the output of HotMAPS.
    It will be used as ground truth for the new clustering method.
    """
    
    # Select only necessary info and rename to match the input file
    cluster = pd.read_csv(out_cluster_path, sep="\t")
    cluster = cluster.copy().rename(columns={"CRAVAT Res" : "Pos"})[["Pos", "Ref AA", "HUGO Symbol", "Min p-value", "q-value"]]
    cluster = cluster.rename(columns={"HUGO Symbol" : "Gene", "Ref AA" : "WT"})
    
    return cluster


## Handle proteins fragments

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
    uniq_pos_barcodes = [len(pos_barcodes[[pos in np.where(cmap[i-1])[0] for pos in pos_barcodes.Pos]].Barcode.explode().unique()) for i in pos_barcodes.Pos]
    pos_barcodes["Samples_in_vol"] = np.array(uniq_pos_barcodes) / tot_samples
    
    ## TO DELETE #####################################################################################################
    # pos_barcodes["Tot_samples"] = tot_samples               
    # pos_barcodes["# barcodes"] = uniq_pos_barcodes
    
    return tot_samples, pos_barcodes


def add_samples_info(mut_gene_df, result_pos_df, samples_in_vol, tot_samples):
    """
    Add information about the ratio of unique samples in the volume of 
    each mutated residues and in each detected community (meta-cluster) 
    to the residues-level output of the tool.
    """

    # Add samples in vol
    result_pos_df = result_pos_df.merge(samples_in_vol.drop(columns=["Barcode"]), on="Pos", how="outer")
    
    # Get per-community ratio of mutated samples
    community_samples = result_pos_df.groupby("Community").apply(lambda x: 
        len(mut_gene_df[[pos in x.Pos.values for pos in mut_gene_df.Pos]].Tumor_Sample_Barcode.unique()))
    community_samples = community_samples / tot_samples
    community_samples = community_samples.reset_index().rename(columns={0 : "Samples_in_comm"})
    
    
    # ## TOP DDELÒETE À#####################################################################################################
    # community_samples["# samples comm"] = result_pos_df.groupby("Community").apply(lambda x: 
    #     len(mut_gene_df[[pos in x.Pos.values for pos in mut_gene_df.Pos]].Tumor_Sample_Barcode.unique()))
    # community_samples["tot # samples comm"] = tot_samples
    
    # Add to residues-level result
    result_pos_df = result_pos_df.merge(community_samples, on="Community", how="outer")
    #result_pos_df.insert(6, "Samples_in_vol", result_pos_df.pop("Samples_in_vol"))
    
    return result_pos_df


def add_nan_clust_cols(result_gene):
    """
    Add columns showing clustering results with only nan for 
    genes that are not tested (not enough mutations, etc).
    """

    result_gene = result_gene.copy()
    
    columns = ["pval", "qval", "C_gene", "C_pos", "C_community",
               "Top_samples_in_vol", "Top_samples_in_comm", 
               "Top_ratio_obs_sim", #"Top_diff_obs_sim", 
               "Clust_mut", "Top_mut_in_vol", "F", 
               "Mut_in_top_F", "Top_F"]
    
    for col in columns:
        result_gene[col] = np.nan

    return result_gene


def sort_cols(result_gene):
    """
    Simply change the order of columns of the genes-level result.
    """

    cols = ['Gene', 'Uniprot_ID', 
            'pval', 'qval', 'C_gene', 'C_pos', 'C_community',
            'Top_samples_in_vol', 'Top_samples_in_comm', 
            'Top_ratio_obs_sim', 'Clust_mut', 'Top_mut_in_vol', 
            'Mut_in_gene', 'F', 'Mut_in_top_F', 'Top_F', 'Status', 
            'Cancer', 'Cohort']

    return result_gene[[col for col in cols if col in result_gene.columns]]