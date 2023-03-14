"""
Contains function to process the experimental p-values.
"""

from scipy.stats import rankdata


def fdr(p_vals):
    """
    Compute false discovery rate.
    """

    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1

    return fdr


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


def get_final_gene_result(result_pos, result_gene, alpha_gene=0.05):
    """
    Output the final dataframe including gene global pval, 
    significant positions, communities, processing status, etc.
    """
    
    # Get significant positions and communities for each gene
    clusters = result_pos[result_pos["C"] == 1].groupby("Gene").apply(lambda x: (x["Pos"].values)).reset_index().rename(
        columns={0 : "C_pos"})
    clusters["C_community"] = result_pos[result_pos["C"] == 1].groupby("Gene").apply(lambda x: x["Community"].values).reset_index(drop=True)
    
    # Annotate each gene with significant hits
    result_gene = clusters.merge(result_gene, on="Gene", how="outer")
    
    # Get gene pval, qval, and largest density omong the most significant hits
    gene_pvals = result_pos.groupby("Gene").apply(lambda x: min(x["pval"].values)).reset_index().rename(
        columns={0 : "pval"})
    gene_pvals["qval"] = fdr(gene_pvals["pval"])
    gene_pvals["Mut_in_vol"] = result_pos.groupby("Gene").apply(lambda x: max(x[x["pval"] == min(x["pval"])].Mut_in_vol)).values

    # Combine pval to significant annotation
    result_gene = gene_pvals.merge(result_gene, on="Gene", how="outer")

    # Add gene binary clustering label and sort genes
    gene_label = result_gene.apply(lambda x: 1 if x.qval < alpha_gene else 0, axis=1)
    result_gene.insert(4, "C_gene", gene_label)
    result_gene.insert(8, "Max_mut_in_vol", result_gene.pop("Mut_in_vol"))
    result_gene.insert(1, "Uniprot_ID", result_gene.pop("Uniprot_ID"))
    result_gene = sort_by_pval_and_mut(result_gene)
    
    return result_gene