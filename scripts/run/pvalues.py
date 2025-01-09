"""
Contains function to process the experimental p-values.
"""

import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests


def fdr(p_vals, alpha=0.05):
    """
    Compute false discovery rate using Benjamini-Hochberg method.
    """

    return multipletests(p_vals, alpha=alpha, method='fdr_bh', is_sorted=True)[1]   


def get_top_vol_info(gene_result_pos):
    """
    Get score, mutations count and other info of the top volume 
    (most significant one) across the mutated ones in each gene. 
    """
    
    if len(gene_result_pos) > 1:
        lowest_pval = gene_result_pos[gene_result_pos['pval'] == gene_result_pos['pval'].min()]
        if len(lowest_pval) > 1:
            top_vol = lowest_pval[lowest_pval['Score_obs_sim'] == lowest_pval['Score_obs_sim'].max()].iloc[0]
        else:
            top_vol = lowest_pval.iloc[0]
    else:
        top_vol = gene_result_pos.iloc[0]
    
    pos_top_vol = top_vol.Pos
    mut_in_top_vol = np.round(top_vol.Mut_in_vol, 2)
    mut_in_top_cl_vol = np.round(top_vol.Mut_in_cl_vol, 2)
    score_obs_sim_top_vol = top_vol.Score_obs_sim
    pae_top_vol = np.round(top_vol.PAE_vol, 2)
    plddt_top_vol = top_vol.pLDDT_vol
    pLDDT_top_cl_vol = np.round(top_vol.pLDDT_cl_vol, 2)
    
    return pos_top_vol, mut_in_top_vol, mut_in_top_cl_vol, score_obs_sim_top_vol, pae_top_vol, plddt_top_vol, pLDDT_top_cl_vol


def get_final_gene_result(result_pos, result_gene, alpha_gene=0.05, sample_info=False):
    """
    Output the final dataframe including gene global pval, qval,
    significant positions, clumps, processing status, etc.
    """

    pos_hits = result_pos[result_pos["C"] == 1]

    if len(pos_hits) > 0:
        # Get significant positions and communities for each gene
        clumps = pos_hits.groupby("Gene").apply(lambda x: (x["Pos"].values)).reset_index().rename(columns={0 : "C_pos"})
        clumps["C_label"] = pos_hits.groupby("Gene").apply(lambda x: x["Clump"].values).reset_index(drop=True)
        # Annotate each gene with significant hits
        result_gene = clumps.merge(result_gene, on="Gene", how="outer")
    else:
        result_gene["C_pos"] = np.nan
        result_gene["C_label"] = np.nan
    
    # Gene pval
    gene_pvals = result_pos.groupby("Gene").apply(lambda x: min(x["pval"].values)).reset_index().rename(columns={0 : "pval"})
 
    # Top volume info
    gene_top_vol_info = result_pos.groupby("Gene").apply(lambda x: get_top_vol_info(x)).apply(pd.Series)
    gene_top_vol_info.columns = ["Pos_top_vol",
                                 "Mut_in_top_vol", 
                                 "Mut_in_top_cl_vol", 
                                 "Score_obs_sim_top_vol", 
                                 "PAE_top_vol", 
                                 "pLDDT_top_vol", 
                                 "pLDDT_top_cl_vol"]
    gene_top_vol_info = gene_top_vol_info.reset_index()
    gene_pvals = gene_pvals.merge(gene_top_vol_info, on="Gene")

    # Sort positions and get qval
    gene_pvals = gene_pvals.sort_values(["pval", "Score_obs_sim_top_vol"], ascending=[True, False]).reset_index(drop=True)
    not_processed_genes_count = sum(~result_gene.Status.str.contains("Processed", na=False))
    gene_pvals["qval"] = fdr(np.concatenate((gene_pvals["pval"], np.repeat(1, not_processed_genes_count))))[:len(gene_pvals)]
                     
    # Combine gene-level clustering result, add label, sort genes, add fragment info
    result_gene = gene_pvals.merge(result_gene, on="Gene", how="outer")
    result_gene["C_gene"] = result_gene.apply(lambda x: 1 if x.qval < alpha_gene else 0, axis=1)
    # result_gene = result_gene.sort_values(["pval", "Score_obs_sim_top_vol"], ascending=[True, False])
    
    # Convert C_pos and C_label to str
    result_gene["C_pos"] = result_gene["C_pos"].apply(lambda x: str(x) if isinstance(x, list) else x)
    result_gene["C_label"] = result_gene["C_label"].apply(lambda x: str(x) if isinstance(x, list) else x)

    return result_gene