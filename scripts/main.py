""" 
The module includes the the main function to run an HotMAPs-inspired 
method that uses the the canonical predicted structure stored in 
AlphaFold db to perform 3D-clustering of mutations using simulations 
and rank based comparison.

###################################### EXAMPLE USAGE ############################################

time python3 main.py -i ../tests/input/HARTWIG_WGS_PANCREAS_2020.in.maf -o ../tests/output/ \
-p ../tests/input/HARTWIG_WGS_PANCREAS_2020.mutrate.json -H 0 -t PANCREAS -C HARTWIG_WGS_PANCREAS_2020

python3 main.py -i ../../../evaluation/datasets/input/maf/PCAWG_WGS_COLORECT_ADENOCA.in.maf \         
-o ../../../evaluation/output \
-p ../../../evaluation/datasets/input/mut_profile/prior.weighted.normalized.PCAWG_WGS_COLORECT_ADENOCA.json \
-s ../datasets/seq_for_mut_prob.csv \
-c ../datasets/cmaps/ \
-u ../datasets/af_uniprot_to_gene_id.json \
-n 10000 \
-H 0 \
-t COREAD \
-C PCAWG_WGS_COLORECT_ADENOCA

#################################################################################################
"""


import json
import numpy as np
import pandas as pd
from progressbar import progressbar
import argparse
from utils.utils import parse_maf_input
from utils.miss_mut_prob import mut_rate_vec_to_dict, get_miss_mut_prob_dict
from utils.clustering import clustering_3d, clustering_3d_frag
from utils.pvalues import get_final_gene_result


def main():
    """
    Wrapper function.
    """

    ## Initialize
    version = "2023_v1.3"    # LAST CHANGE: obtain HUGO mapping from seq_df

    # Parser
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)

    parser.add_argument("-i", "--input_maf", help="Path of the maf file used as input", type=str, required=True)
    parser.add_argument("-o", "--output_dir", help="Path to output directory", type=str, required=True)

    group.add_argument("-p", "--mut_profile", help="Path to the mut profile (list of 96 floats) of the cohort (json)", type=str)
    group.add_argument("-P", "--miss_mut_prob", help="Path to the dict of missense mut prob of each protein based on mut profile of the cohort (json)",  
                       type=str)
    
    parser.add_argument("-s", "--seq_df", 
                        help="Path to the dataframe including DNA and protein seq of all gene/proteins (all AF predicted ones)", 
                        type=str, 
                        default="../datasets/seq_for_mut_prob.csv")       
    parser.add_argument("-c", "--cmap_path", help="Path to the directory containting the contact map of each protein", type=str, 
                        default="../datasets/cmaps/")

    parser.add_argument("-n", "--n_iterations", help="Number of densities to be simulated", type=int, default=10000)
    parser.add_argument("-a", "--alpha_level_res", help="Significant threshold for the p-value of protein residues", type=float, default=0.01)
    parser.add_argument("-A", "--alpha_level_gene", help="Significant threshold for the global q-value of the genes", type=float, default=0.05)
    parser.add_argument("-H", "--hits_only", help="if 1 returns only positions in clusters, if 0 returns all", type=int, default=1)
    
    parser.add_argument("-t", "--cancer_type", help="Cancer type", type=str)
    parser.add_argument("-C", "--cohort_name", help="Name of the cohort", type=str)
    args = parser.parse_args()

    maf_input_path = args.input_maf
    mut_profile_path = args.mut_profile
    miss_mut_prob_path = args.miss_mut_prob
    seq_df_path = args.seq_df
    cmap_path = args.cmap_path
    num_iteration = args.n_iterations
    output_dir = args.output_dir
    cohort = args.cohort_name
    cancer_type = args.cancer_type
    alpha_pos = args.alpha_level_res
    alpha_gene = args.alpha_level_gene
    hits_only = args.hits_only

    print(f"Starting 3D-clustering [{version}]..\n")
    print(f"Iterations: {num_iteration}")
    print(f"Significant level (position): {alpha_pos}")
    print(f"Global significant level (gene): {alpha_gene}")
    print(f"Cohorts: {cohort}")
    print(f"Cancer type: {cancer_type}")
    print(f"Output directory: {output_dir}")


    ## Load input and df of DNA sequences

    # MAF input
    data = parse_maf_input(maf_input_path)

    # Seq df for missense mut prob
    seq_df = pd.read_csv(seq_df_path)


    ## Run

    result_pos_lst = []
    result_gene_lst = []

    # Get genes with enough mut
    genes = data.groupby("Gene").apply(lambda x: len(x))
    genes_mut = genes[genes >= 2].index
    genes_no_mut = genes[genes < 2].index
    result_gene = pd.DataFrame({"Gene" : genes_no_mut,
                                "Uniprot_ID" : np.nan,
                                "No_frag" : np.nan,
                                "Mut_in_gene" : 1,
                                "Max_mut_pos" : np.nan,
                                "Structure_max_pos" : np.nan,
                                "Status" : "No_mut"})
    result_gene_lst.append(result_gene)

    # Get genes with corresponding Uniprot-ID mapping
    gene_to_uniprot_dict = {gene : uni_id for gene, uni_id in seq_df[["Gene", "Uniprot_ID"]].drop_duplicates().values}
    genes_mapped = [gene for gene in genes_mut if gene in gene_to_uniprot_dict.keys()]
    seq_df = seq_df[[gene in genes_mapped for gene in seq_df["Gene"]]].reset_index(drop=True)
    genes_no_mapping = genes[[gene in genes_mut and gene not in gene_to_uniprot_dict.keys() for gene in genes.index]]
    result_gene = pd.DataFrame({"Gene" : genes_no_mapping.index,
                                "Uniprot_ID" : np.nan,
                                "No_frag" : np.nan,
                                "Mut_in_gene" : genes_no_mapping.values,
                                "Max_mut_pos" : np.nan,
                                "Structure_max_pos" : np.nan,
                                "Status" : "No_ID_mapping"})
    result_gene_lst.append(result_gene)

    # Missense mut prob  
    if miss_mut_prob_path is not None:
        # Load dict with miss prob of each prot
        miss_prob_dict = json.load(open(miss_mut_prob_path))
    else:
        # Compute dict from mut profile of the cohort and dna sequences
        mut_profile = json.load(open(mut_profile_path))
        print(f"\nComputing missense mut probabilities, # proteins/fragment: {len(seq_df)}")
        mut_profile = mut_rate_vec_to_dict(mut_profile)
        miss_prob_dict = get_miss_mut_prob_dict(mut_rate_dict=mut_profile, seq_df=seq_df)

    # Process gene
    print("Performing 3D clustering on genes with enough mutations..")
    for gene in progressbar(genes_mapped):
        mut_gene_df = data[data["Gene"] == gene]
        uniprot_id = gene_to_uniprot_dict[gene]   

        # If there is a single fragment
        if seq_df[seq_df["Gene"] == gene].F.max() == 1:      
            
            try:
                pos_result, result_gene = clustering_3d(gene,
                                                         uniprot_id, 
                                                         mut_gene_df, 
                                                         cmap_path,
                                                         miss_prob_dict,
                                                         fragment=1,
                                                         alpha=alpha_pos,
                                                         num_iteration=num_iteration,
                                                         hits_only=hits_only)
                result_gene_lst.append(result_gene)
                if pos_result is not None:
                    result_pos_lst.append(pos_result)
                                                              
            except:                                                     # >>>> Should raise a better exception to capture a more specific error
                result_gene = pd.DataFrame({"Gene" : gene,
                                            "Uniprot_ID" : uniprot_id,
                                            "No_frag" : np.nan,
                                            "Mut_in_gene" : np.nan,
                                            "Max_mut_pos" : np.nan,
                                            "Structure_max_pos" : np.nan,
                                            "Status" : "Not_processed"},
                                            index=[0])
                result_gene_lst.append(result_gene)

        # If the protein is fragmented
        else:
            
            try:
                pos_result, result_gene = clustering_3d_frag(gene, 
                                                            uniprot_id,
                                                            mut_gene_df,
                                                            cmap_path,
                                                            miss_prob_dict,
                                                            alpha=alpha_pos,
                                                            num_iteration=num_iteration,
                                                            hits_only=hits_only)
                result_gene_lst.append(result_gene)
                if pos_result is not None:
                    result_pos_lst.append(pos_result)

            except:
                result_gene = pd.DataFrame({"Gene" : gene,
                                            "Uniprot_ID" : uniprot_id,
                                            "No_frag" : np.nan,
                                            "Mut_in_gene" : np.nan,
                                            "Max_mut_pos" : np.nan,
                                            "Structure_max_pos" : np.nan,
                                            "Status" : "Not_processed_frag"},
                                            index=[0])
                result_gene_lst.append(result_gene)


    ## Save 
    print("\nSaving..")
    result_gene = pd.concat(result_gene_lst)
    result_gene["Cancer"] = cancer_type
    result_gene["Cohort"] = cohort

    if len(result_pos_lst) == 0:
        print(f"Did not processed any genes\n")
        result_gene.to_csv(f"{output_dir}/{cancer_type}.{cohort}.3d_clustering_genes.csv", index=False)
    else:
        result_pos = pd.concat(result_pos_lst)
        result_pos["Cancer"] = cancer_type
        result_pos["Cohort"] = cohort
        result_pos.to_csv(f"{output_dir}/{cancer_type}.{cohort}.3d_clustering_pos.csv", index=False)

        # Get gene global pval, qval, and clustering annotations
        result_gene = get_final_gene_result(result_pos, result_gene, alpha_gene)
        result_gene = result_gene.drop(columns = ["Max_mut_pos", "Structure_max_pos"]) # comment out if willing to look at isoforms incongruence
        result_gene.to_csv(f"{output_dir}/{cancer_type}.{cohort}.3d_clustering_genes.csv", index=False)


if __name__ == "__main__":
    main()
 