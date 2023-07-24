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


def get_final_gene_result(result_pos, result_gene, alpha_gene=0.05):
    """
    Output the final dataframe including gene global pval, s
    significant positions, clusters, processing status, etc.
    """

    pos_hits = result_pos[result_pos["C"] == 1]

    if len(pos_hits) > 0:

        # Get significant positions and communities for each gene
        clusters = pos_hits.groupby("Gene").apply(lambda x: (x["Pos"].values)).reset_index().rename(columns={0 : "C_pos"})
        clusters["C_label"] = pos_hits.groupby("Gene").apply(lambda x: x["Cluster"].values).reset_index(drop=True)
        
        # Annotate each gene with significant hits
        result_gene = clusters.merge(result_gene, on="Gene", how="outer")

    else:
        result_gene["C_pos"] = np.nan
        result_gene["C_label"] = np.nan
    
    # Gene pval, qval
    gene_pvals = result_pos.groupby("Gene").apply(lambda x: min(x["pval"].values)).reset_index().rename(columns={0 : "pval"})
    
    # Sample info and largest density among hits
    # > NB: samples info for fragments will be displayed as they are individual proteins <
    gene_pvals["Tot_samples"] = result_pos.groupby("Gene").apply(lambda x: x["Tot_samples"].unique()[0]).values
    gene_pvals["Samples_in_top_vol"] = result_pos.groupby("Gene").apply(lambda x: max(x[x["pval"] == min(x["pval"])].Samples_in_vol)).values
    gene_pvals["Samples_in_top_cl_vol"] = result_pos.groupby("Gene").apply(lambda x: max(x[x["pval"] == min(x["pval"])].Samples_in_cl_vol)).values
    gene_pvals["Mut_in_top_cl_vol"] = result_pos.groupby("Gene").apply(lambda x: max(x[x["pval"] == min(x["pval"])].Mut_in_cl_vol)).values
    gene_pvals["Ratio_obs_sim_top_vol"] = result_pos.groupby("Gene").apply(lambda x: max(x[x["pval"] == min(x["pval"])].Ratio_obs_sim)).values
    gene_pvals["Mut_in_top_vol"] = result_pos.groupby("Gene").apply(lambda x: max(x[x["pval"] == min(x["pval"])].Mut_in_vol)).values
    gene_pvals["PAE_top_vol"] = result_pos.groupby("Gene").apply(lambda x: max(x[x["pval"] == min(x["pval"])].PAE_vol)).values
    gene_pvals["pLDDT_top_vol"] = result_pos.groupby("Gene").apply(lambda x: max(x[x["pval"] == min(x["pval"])].pLDDT_vol)).values
    gene_pvals["pLDDT_top_cl_vol"] = result_pos.groupby("Gene").apply(lambda x: max(x[x["pval"] == min(x["pval"])].pLDDT_cl_vol)).values
    
    # Sort positions and get qval
    gene_pvals = gene_pvals.sort_values(["pval"], ascending=True).reset_index(drop=True)
    not_processed_genes_count = sum(~result_gene.Status.str.contains("Processed", na=False))
    gene_pvals["qval"] = fdr(np.concatenate((gene_pvals["pval"], np.repeat(1, not_processed_genes_count))))[:len(gene_pvals)]
                     
    # Combine gene-level clustering result, add label, sort genes, add fragment info
    result_gene = gene_pvals.merge(result_gene, on="Gene", how="outer")
    result_gene["C_gene"] = result_gene.apply(lambda x: 1 if x.qval < alpha_gene else 0, axis=1)
    result_gene = result_gene.sort_values(["pval", "Ratio_obs_sim_top_vol"], ascending=[True, False])

    return result_gene