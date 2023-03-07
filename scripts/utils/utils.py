import re
import pandas as pd
import numpy as np
import csv
import requests
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

## Parsers

def parse_maf_input(maf_input_path):
    """
    Parse in.maf file which is used as 
    input for the HotMAPS method.
    """

    # Load
    maf = pd.read_csv(maf_input_path, sep="\t")

    # Select only missense mutation and extract Gene name and mut
    maf = maf.loc[maf.Variant_Classification == "Missense_Mutation"].copy()
    maf["Pos"] = maf.loc[:, "HGVSp_Short"].apply(lambda x: int(re.sub("\\D", "", (x[2:]))))
    maf["WT"] = maf["HGVSp_Short"].apply(lambda x: re.findall("\\D", x[2:])[0])
    maf["Mut"] = maf["HGVSp_Short"].apply(lambda x: re.findall("\\D", x[2:])[1])
    maf = maf[["Hugo_Symbol", "Pos", "WT", "Mut"]]
    maf = maf.sort_values("Pos").rename(columns={"Hugo_Symbol" : "Gene"})
    
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


## Eval    >>>> DEPRECATED <<<<< 

def get_evaluation_df(test_result, mut_gene_df, drop_extra=False):
    """
    Merge the results from the statistical test (excess of densities 
    of mutations) with the score assigned from HotMAPs.
    """

    if drop_extra:
        cols_to_drop = ["Gene", "WT", "Mut"]
        mut_gene_df = mut_gene_df.drop(columns = cols_to_drop)
    
    # Merge the new result (mutations predicted to be in positions enriched of mutations)
    # with the original data (mutations that were not predicted to be enriched will have NA in that col)
    mut_gene_df = mut_gene_df.sort_values("Pos")
    evaluation_df = test_result.merge(mut_gene_df, how = "outer").sort_values("Pos").reset_index(drop=True)
    
    return evaluation_df


def evaluate(evaluation_df):

    evaluation_df = evaluation_df[["Gene", 
                                   "Pos", 
                                   "New_hotmaps", 
                                   "HotMAPs"]].drop_duplicates(keep="first")
    
    all_mut = len(evaluation_df)

    all_pos = evaluation_df[evaluation_df["HotMAPs"] == 1]
    pos_pred = evaluation_df[evaluation_df["New_hotmaps"] == 1]
    tp_pred = pos_pred[pos_pred["HotMAPs"] == 1]
    fp_pred = pos_pred[pos_pred["HotMAPs"] == 0]

    all_neg = evaluation_df[evaluation_df["HotMAPs"] == 0]
    neg_pred = evaluation_df[evaluation_df["New_hotmaps"] == 0]
    tn_pred = neg_pred[neg_pred["HotMAPs"] == 0]
    fn_pred = neg_pred[neg_pred["HotMAPs"] == 1]
    
    if len(all_pos) > 0:
        tpr = round(len(tp_pred) / (len(tp_pred) + len(fn_pred)), 3)
        if len(pos_pred) > 0:
            precision = round(len(tp_pred) / (len(tp_pred) + len(fp_pred)), 3)
        else:
            precision = np.nan
    else:
        tpr = np.nan
        precision = np.nan
    if len(all_neg):
        tnr = round(len(tn_pred) / (len(tn_pred) + len(fp_pred)), 3)
    else:
        tnr = np.nan

    return pd.DataFrame({"Gene" : evaluation_df.Gene.unique()[0],
                         "All_mut" : all_mut, 
                         "All_pos" : len(all_pos),
                         "Pos_pred" : len(pos_pred), 
                         "TP_pred" : len(tp_pred), 
                         "FP_pred" : len(fp_pred), 
                         "All_neg" : len(all_neg),
                         "Neg_pred" : len(neg_pred),
                         "TN_pred" : len(tn_pred),
                         "FN_pred" : len(fn_pred),
                         "TPR" : tpr, 
                         "Precis" : precision,
                         "TNR" : tnr,}, 
                        index=[1])

    
def plot_venn(eval_df_final, save=False, path="venn.png"):
    """
    Venn diagram for genes found by HotMAPS, 
    the New HotMAPS, and both methods.
    """
    
    both = len(eval_df_final[(eval_df_final.New_hotmaps == 1) & (eval_df_final.HotMAPs == 1)].Gene.unique())
    new_hotmaps = len(eval_df_final[eval_df_final.New_hotmaps == 1].Gene.unique()) - both
    hotmaps = len(eval_df_final[eval_df_final.HotMAPs == 1].Gene.unique()) - both

    venn2(subsets = (hotmaps, new_hotmaps, both), set_labels = ('HotMAPS', 'New HotMAPS'))
    plt.title("Genes with one or more clusters", fontsize=14)
    if save: 
        plt.savefig(path, bbox_inches="tight", dpi=300)
        plt.clf()
    else:
        plt.show()