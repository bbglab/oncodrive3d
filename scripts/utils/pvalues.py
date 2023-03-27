"""
Contains function to process the experimental p-values.
"""

from statsmodels.stats.multitest import multipletests
import numpy as np


def fdr(p_vals, alpha=0.05):
    """
    Compute false discovery rate using Benjamini-Hochberg method.
    """

    return multipletests(p_vals, alpha=alpha, method='fdr_bh', is_sorted=True)[1]


def add_nan_clust_cols(result_gene):
    """
    Add columns showing clustering results with only nan for 
    genes that are not tested (not enough mutations, etc).
    """

    result_gene = result_gene.copy()
    result_gene.insert(2, "pval", np.nan)
    result_gene.insert(3, "qval", np.nan)
    result_gene.insert(4, "C_gene", np.nan)
    result_gene.insert(5, "C_pos", np.nan)
    result_gene.insert(6, "C_community", np.nan)
    result_gene.insert(7, "Top_ratio_obs_sim", np.nan)
    #result_gene.insert(8, "Top_diff_obs_sim", np.nan)
    result_gene.insert(8, "Clust_mut", np.nan)
    result_gene.insert(8, "Top_mut_in_vol", np.nan)
    result_gene.insert(11, "F", result_gene.pop("F"))
    result_gene.insert(12, "Mut_in_top_F", np.nan)
    result_gene.insert(13, "Top_F", np.nan)

    return result_gene


def add_gene_binary_and_sort(result_gene, alpha_gene):
    """
    Add gene binary clustering label, re-order columns for the 
    final gene-level, and sort genes according to p-values and 
    the ratio between observed and simulated score.
    """

    # Assign binary label
    gene_label = result_gene.apply(lambda x: 1 if x.qval < alpha_gene else 0, axis=1)

    # Re-order columns
    result_gene.insert(1, "Uniprot_ID", result_gene.pop("Uniprot_ID"))
    result_gene.insert(9, "Top_ratio_obs_sim", result_gene.pop("Top_ratio_obs_sim"))
    #result_gene.insert(9, "Top_diff_obs_sim", result_gene.pop("Top_diff_obs_sim"))
    result_gene.insert(10, "Clust_mut", result_gene.pop("Clust_mut"))
    result_gene.insert(10, "Top_mut_in_vol", result_gene.pop("Top_mut_in_vol"))
    result_gene.insert(10, "Mut_in_gene", result_gene.pop("Mut_in_gene"))   
    result_gene.insert(10, "F", result_gene.pop("F"))   
    result_gene.insert(4, "C_gene", gene_label)     

    # Sort genes
    result_gene = result_gene.sort_values(['pval', 'Top_ratio_obs_sim'], ascending=[True, False])

    return result_gene


def add_frag_info(result_pos, result_gene):
    """
    Add total mutation in the fragment and number fragment with the most significant cluster.
    """

    frag_info = result_pos.groupby("Gene").apply(lambda x: 
                                                 max(x[x["pval"] == min(x["pval"])].Mut_in_gene) if max(x.F) > 1 
                                                 else np.nan).reset_index().rename(columns = {0 : "Mut_in_top_F"}) 
    frag_info["Top_F"] = result_pos.groupby("Gene").apply(lambda x: max(x[x["pval"] == min(x["pval"])].F) if max(x.F) > 1 else np.nan).values                                                  
    result_gene = result_gene.merge(frag_info, on = "Gene", how = "outer")
    result_gene.insert(13, "Mut_in_top_F", result_gene.pop("Mut_in_top_F"))
    result_gene.insert(14, "Top_F", result_gene.pop("Top_F"))

    return result_gene


def get_final_gene_result(result_pos, result_gene, alpha_gene=0.05):
    """
    Output the final dataframe including gene global pval, 
    significant positions, communities, processing status, etc.
    """

    pos_hits = result_pos[result_pos["C"] == 1]

    if len(pos_hits) > 0:

        # Get significant positions and communities for each gene
        clusters = pos_hits.groupby("Gene").apply(lambda x: (x["Pos"].values)).reset_index().rename(columns={0 : "C_pos"})
        clusters["C_community"] = pos_hits.groupby("Gene").apply(lambda x: x["Community"].values).reset_index(drop=True)
        
        # Annotate each gene with significant hits
        result_gene = clusters.merge(result_gene, on="Gene", how="outer")

    else:
        result_gene.insert(0, "C_pos", np.nan)
        result_gene.insert(1, "C_community", np.nan)
    
    # Get gene pval, qval, and largest density among the most significant hits
    gene_pvals = result_pos.groupby("Gene").apply(lambda x: min(x["pval"].values)).reset_index().rename(columns={0 : "pval"})
    gene_pvals["Top_ratio_obs_sim"] = result_pos.groupby("Gene").apply(lambda x: max(x[x["pval"] == min(x["pval"])].Ratio_obs_sim)).values
    #gene_pvals["Top_diff_obs_sim"] = result_pos.groupby("Gene").apply(lambda x: max(x[x["pval"] == min(x["pval"])].Diff_obs_sim)).values
    gene_pvals["Top_mut_in_vol"] = result_pos.groupby("Gene").apply(lambda x: max(x[x["pval"] == min(x["pval"])].Mut_in_vol)).values

    gene_pvals = gene_pvals.sort_values(["pval"], ascending=True).reset_index(drop=True)
    gene_pvals["qval"] = fdr(gene_pvals["pval"])

    # Combine gene-level clustering result, add label, sort genes, add fragment info
    result_gene = gene_pvals.merge(result_gene, on="Gene", how="outer")
    result_gene = add_gene_binary_and_sort(result_gene, alpha_gene)
    result_gene = add_frag_info(result_pos, result_gene)

    return result_gene