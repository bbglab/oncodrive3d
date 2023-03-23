"""
Contains function to process the experimental p-values.
"""

from statsmodels.stats.multitest import multipletests
import numpy as np
#from scipy.stats import rankdata


def fdr(p_vals, alpha=0.05):
    """
    Compute false discovery rate using Benjamini-Hochberg method.
    """

    return multipletests(p_vals, alpha=alpha, method='fdr_bh', is_sorted=True)[1]


# def fdr(p_vals, alpha=0.05):
#     """
#     Compute false discovery rate.
#     """

#     ranked_p_values = rankdata(p_vals)
#     fdr = p_vals * len(p_vals) / ranked_p_values
#     fdr[fdr > 1] = 1

#     return fdr


def sort_by_pval_and_mut(result_gene):
    """
    Sort by p-value and break the tie using the 
    density of mutations (larger the better) and 
    total mutation in the gene (smaller the better). 
    """
    
    result_gene = result_gene.copy()
    result_gene["Max_mut_in_vol"] = -result_gene["Max_mut_in_vol"]
    result_gene = result_gene.sort_values(["pval", "Max_mut_in_vol", "Mut_in_gene"], ascending=True).reset_index(drop=True)
    result_gene["Max_mut_in_vol"] = -result_gene["Max_mut_in_vol"]
    
    return result_gene


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
    result_gene.insert(6, "Max_mut_in_vol", np.nan)

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
    gene_pvals["Max_mut_in_vol"] = result_pos.groupby("Gene").apply(lambda x: max(x[x["pval"] == min(x["pval"])].Mut_in_vol)).values
    gene_pvals = gene_pvals.sort_values(["pval"], ascending=True).reset_index(drop=True)
    gene_pvals["qval"] = fdr(gene_pvals["pval"])

    # Combine pval to significant annotation
    result_gene = gene_pvals.merge(result_gene, on="Gene", how="outer")

    # Add gene binary clustering label and sort genes
    gene_label = result_gene.apply(lambda x: 1 if x.qval < alpha_gene else 0, axis=1)
    result_gene.insert(4, "C_gene", gene_label)
    result_gene.insert(8, "Max_mut_in_vol", result_gene.pop("Max_mut_in_vol"))
    result_gene.insert(1, "Uniprot_ID", result_gene.pop("Uniprot_ID"))
    result_gene = sort_by_pval_and_mut(result_gene)

    return result_gene