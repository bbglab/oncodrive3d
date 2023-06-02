#!/usr/bin/env python

""" 
The module includes the the main function to run an HotMAPs-inspired 
method that uses the the canonical predicted structure stored in 
AlphaFold db to perform 3D-clustering of mutations using simulations 
and rank based comparison.

###################################### EXAMPLE USAGE ############################################

time python3 main.py -i ../tests/input/HARTWIG_WGS_PANCREAS_2020.in.maf -o ../tests/output/ \
-p ../tests/input/HARTWIG_WGS_PANCREAS_2020.mutrate.json -H 0 -t PANCREAS -C HARTWIG_WGS_PANCREAS_2020 -e 1

time python3 main.py -i /workspace/projects/clustering_3d/evaluation/datasets/input/maf/ICGC_WXS_HCC_LICA_CN_STRELKA_2019.in.maf \
-o /workspace/projects/clustering_3d/dev_testing/output/ \
-p /workspace/projects/clustering_3d/evaluation/datasets/input/mut_profile/ICGC_WXS_HCC_LICA_CN_STRELKA_2019.mutrate.json \
-H 0 -t PANCREAS -C ICGC_WXS_HCC_LICA_CN_STRELKA_2019

time python3 main.py -i /workspace/projects/clustering_3d/evaluation/datasets/input/maf/STJUDE_WGS_D_LGG_2018.in.maf \
-o /workspace/projects/clustering_3d/dev_testing/output -p /workspace/projects/clustering_3d/evaluation/datasets/input/mut_profile/STJUDE_WGS_D_LGG_2018.mutrate.json \
-H 0 -t LGG -C STJUDE_WGS_D_LGG_2018


time python3 main.py -i /workspace/projects/clustering_3d/evaluation/datasets/input/maf/ICGC_WXS_AML_LAML_KR_VARSCAN_2019.in.maf \
-o /workspace/projects/clustering_3d/dev_testing/output/ \
-p /workspace/projects/clustering_3d/evaluation/datasets/input/mut_profile/ICGC_WXS_AML_LAML_KR_VARSCAN_2019.mutrate.json \
-H 0 -C ICGC_WXS_AML_LAML_KR_VARSCAN_2019

time python3 main.py -i /workspace/projects/clustering_3d/evaluation/datasets/input/maf/STJUDE_WGS_D_ALL_2018.in.maf \
-o /workspace/projects/clustering_3d/dev_testing/output/ \
-p /workspace/projects/clustering_3d/evaluation/datasets/input/mut_profile/STJUDE_WGS_D_ALL_2018.mutrate.json \
-H 0 -C STJUDE_WGS_D_ALL_2018

time python3 main.py -i /workspace/projects/clustering_3d/evaluation/datasets/input/maf/HARTWIG_WGS_PLMESO_2020.in.maf \
-o /workspace/projects/clustering_3d/dev_testing/output/ \
-p /workspace/projects/clustering_3d/evaluation/datasets/input/mut_profile/HARTWIG_WGS_PLMESO_2020.mutrate.json \
-H 0 -C HARTWIG_WGS_PLMESO_2020 -e 1

time python3 main.py -i /workspace/projects/clustering_3d/evaluation/datasets/input/maf/HARTWIG_WGS_NSCLC_2020.in.maf \
-o /workspace/projects/clustering_3d/dev_testing/output/ \
-p /workspace/projects/clustering_3d/evaluation/datasets/input/mut_profile/HARTWIG_WGS_NSCLC_2020.mutrate.json \
-H 0 -C HARTWIG_WGS_NSCLC_2020 -e 1

## CH

time python3 main.py -i /workspace/projects/clustering_3d/evaluation/datasets/datasets_ch/input/maf/OTHER_WXS_CH_IMPACT_PANEL.in.maf \
-o /workspace/projects/clustering_3d/dev_testing/output/ \
-p /workspace/projects/clustering_3d/evaluation/datasets/datasets_ch/input/mut_profile/OTHER_WXS_CH_IMPACT_PANEL.mutrate.json \
-H 0 -t CH -C OTHER_WXS_CH_IMPACT_PANEL -e 1

time python3 main.py -i /workspace/projects/clustering_3d/evaluation/datasets/datasets_ch/input/maf/OTHER_WXS_TCGA_FULL.in.maf \
-o /workspace/projects/clustering_3d/evaluation/tool_output/run_20230512_ch/ \
-p /workspace/projects/clustering_3d/evaluation/datasets/datasets_ch/input/mut_profile/OTHER_WXS_TCGA_FULL.mutrate.json \
-H 0 -t CH -C OTHER_WXS_TCGA_FULL -e 1

#################################################################################################
"""


import argparse
import json
import numpy as np
import pandas as pd
import os
from datetime import datetime
from progressbar import progressbar
from utils.utils import parse_maf_input, add_nan_clust_cols, sort_cols
from utils.miss_mut_prob import mut_rate_vec_to_dict, get_miss_mut_prob_dict
from utils.clustering import clustering_3d, clustering_3d_frag
from utils.pvalues import get_final_gene_result


def init_parser():
    """
    Initialize parser for the main function.
    """

    ## Parser
    parser = argparse.ArgumentParser()

    # Required input
    parser.add_argument("-i", "--input_maf", help = "Path of the maf file used as input", type=str, required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-p", "--mut_profile", help = "Path to the mut profile (list of 96 floats) of the cohort (json)", type=str)
    group.add_argument("-P", "--miss_mut_prob", help = "Path to the dict of missense mut prob of each protein based on mut profile of the cohort (json)",  type=str)   ### <<< Probably unnecessary >>>
    
    # Required output
    parser.add_argument("-o", "--output_dir", help = "Path to output directory", type=str, required=True)
    
    # Precomputed files input
    parser.add_argument("-s", "--seq_df", help = "Path to the dataframe including DNA and protein seq of all gene/proteins (all AF predicted ones)", type=str)       
    parser.add_argument("-c", "--cmap_path", help = "Path to the directory containting the contact map of each protein", type=str)
    parser.add_argument("-d", "--plddt_path", help = "Path to the pandas dataframe including the AF model confidence of all proteins", type=str)

    # Parameters
    parser.add_argument("-n", "--n_iterations", help = "Number of densities to be simulated", type=int, default=10000)
    parser.add_argument("-a", "--alpha_level", help = "Significant threshold for the p-value of res and gene", type=float, default=0.01)
    parser.add_argument("-H", "--hits_only", help = "If 1 returns only positions in clusters, if 0 returns all", type=int, default=0)
    parser.add_argument("-e", "--ext_hits", 
                        help = "If 1 extend clusters to all mutated residues in the significant volumes, if 0 extend only to the ones having an anomaly > expected", 
                        type=int, default=1)
    parser.add_argument("-f", "--fragments", help = "Enable processing of fragmented proteins (AF-F)", type=int, default=0)
    
    # Metadata annotation
    parser.add_argument("-t", "--cancer_type", help = "Cancer type", type=str)
    parser.add_argument("-C", "--cohort_name", help = "Name of the cohort", type=str)

    return parser.parse_args()


def main():
    """
    Wrapper function.
    """

    ## Initialize
    version = "v_2023_06_02"    # LAST CHANGE: Enable Mus musculus processing and other organisms

    # Parser
    args = init_parser()

    maf_input_path = args.input_maf
    mut_profile_path = args.mut_profile
    miss_mut_prob_path = args.miss_mut_prob
    seq_df_path = args.seq_df
    cmap_path = args.cmap_path
    num_iteration = args.n_iterations
    output_dir = args.output_dir
    cohort = args.cohort_name
    cancer_type = args.cancer_type
    alpha = args.alpha_level
    hits_only = args.hits_only
    fragments = args.fragments
    ext_hits = args.ext_hits
    plddt_path = args.plddt_path

    dir_path = os.path.abspath(os.path.dirname(__file__))
    if plddt_path is None:
        plddt_path = "/workspace/projects/alphafold_features/feature_extraction/model_confidence/confidence_df_all_prot.csv"
    if cmap_path is None:
        cmap_path = f"{dir_path}/../datasets/cmaps/"
    if seq_df_path is None:
        seq_df_path = f"{dir_path}/../datasets/seq_for_mut_prob.csv"
    if cancer_type is None:
        cancer_type = np.nan
    if cohort is None:
        date = datetime.now()
        date.strftime("%m-%d-%Y_%H-%M-%S")
        cohort = f"cohort_{date}"

    print(f"Starting 3D-clustering [{version}]..\n")
    print(f"Path to contact maps: {cmap_path}")
    print(f"Path to DNA sequences: {seq_df_path}")
    print(f"Path to pLDDT scores: {plddt_path}")
    print(f"Iterations: {num_iteration}")
    print(f"Significant level: {alpha}")
    print(f"Extend hits: {ext_hits}")
    print(f"Fragments: {fragments}")
    print(f"Cohort: {cohort}")
    print(f"Cancer type: {cancer_type}")
    print(f"Output directory: {output_dir}")


    ## Load input and df of DNA sequences

    # MAF input
    data = parse_maf_input(maf_input_path, keep_samples_id=True)
    #data = data[data["Gene"] == "TP53"] #################################################################
    
    # data = data[[g in ("TP53", "CDKN2A") for g in data.Gene]]
    # data = data[[g in ("TP53", "ZNF615", "ZNF816", "COL6A3", "TTN", "ARMC4", "C10orf71", "ATP5MF-PTCD1", "NRAF", "BRAF",
    #                    'ABCA2', 'ABCC11', 'ACCS', 'BAP1', 'BTBD3', "ATM", "APOB", "FLNB", "PNMA8A", "AANAT") for g in data.Gene]]

    # Seq df for missense mut prob
    seq_df = pd.read_csv(seq_df_path)
    
    # Model confidence
    plddt_df = pd.read_csv(plddt_path)


    ## Run

    result_pos_lst = []
    result_gene_lst = []

    # Get genes with enough mut
    genes = data.groupby("Gene").apply(lambda x: len(x))
    genes_mut = genes[genes >= 2]
    genes_no_mut = genes[genes < 2].index
    result_gene = pd.DataFrame({"Gene" : genes_no_mut,
                                "Uniprot_ID" : np.nan,
                                "F" : np.nan,
                                "Mut_in_gene" : 1,
                                "Max_mut_pos" : np.nan,
                                "Structure_max_pos" : np.nan,
                                "Status" : "No_mut"})
    result_gene_lst.append(result_gene)   

    # Get genes with corresponding Uniprot-ID mapping
    gene_to_uniprot_dict = {gene : uni_id for gene, uni_id in seq_df[["Gene", "Uniprot_ID"]].drop_duplicates().values}
    genes_to_process = [gene for gene in genes_mut.index if gene in gene_to_uniprot_dict.keys()]
    seq_df = seq_df[[gene in genes_to_process for gene in seq_df["Gene"]]].reset_index(drop=True)
    genes_no_mapping = genes[[gene in genes_mut.index and gene not in gene_to_uniprot_dict.keys() for gene in genes.index]]
    result_gene = pd.DataFrame({"Gene" : genes_no_mapping.index,
                                "Uniprot_ID" : np.nan,
                                "F" : np.nan,
                                "Mut_in_gene" : genes_no_mapping.values,
                                "Max_mut_pos" : np.nan,
                                "Structure_max_pos" : np.nan,
                                "Status" : "No_ID_mapping"})
    result_gene_lst.append(result_gene)
    
    # Filter on fragmented (AF-F) genes
    if fragments == False:
        # Return the fragmented genes as non processed output
        genes_frag = seq_df.groupby("Gene").F.max()
        genes_frag = genes_frag[genes_frag > 1].index 
        genes_frag_mut = genes_mut[[gene in genes_frag for gene in genes_mut.index]]
        genes_frag_id = [gene_to_uniprot_dict[gene] for gene in genes_frag_mut.index]                  
        result_gene = pd.DataFrame({"Gene" : genes_frag,
                                    "Uniprot_ID" : genes_frag_id,
                                    "F" : np.nan,
                                    "Mut_in_gene" : genes_frag_mut.values,
                                    "Max_mut_pos" : np.nan,
                                    "Structure_max_pos" : np.nan,
                                    "Status" : "Fragmented"})
        result_gene_lst.append(result_gene)
        # Filter out from genes to process and seq df
        genes_to_process = [gene for gene in genes_to_process if gene not in genes_frag]
        seq_df = seq_df[[gene in genes_to_process for gene in seq_df["Gene"]]].reset_index(drop=True)

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
    print("Performing 3D clustering..")
    for gene in progressbar(genes_to_process):
        mut_gene_df = data[data["Gene"] == gene]
        uniprot_id = gene_to_uniprot_dict[gene]
        
        # Add confidence to mut_gene_df
        plddt_df_gene_df = plddt_df[plddt_df["Uniprot_ID"] == uniprot_id]
        mut_gene_df = mut_gene_df.merge(plddt_df_gene_df, on = ["Pos"], how = "left")

        # If there is a single fragment
        if seq_df[seq_df["Gene"] == gene].F.max() == 1:   
            
            try:
                pos_result, result_gene = clustering_3d(gene,
                                                        uniprot_id, 
                                                        mut_gene_df, 
                                                        cmap_path,
                                                        miss_prob_dict,
                                                        fragment=1,
                                                        alpha=alpha,
                                                        num_iteration=num_iteration,
                                                        hits_only=hits_only,
                                                        ext_hits=ext_hits)
                result_gene_lst.append(result_gene)
                if pos_result is not None:
                    result_pos_lst.append(pos_result)
                                                              
            except:                                                     # >>>> Should raise a better exception to capture a more specific error                
                
                result_gene = pd.DataFrame({"Gene" : gene,
                                            "Uniprot_ID" : uniprot_id,
                                            "F" : np.nan,
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
                                                            alpha=alpha,
                                                            num_iteration=num_iteration,
                                                            hits_only=hits_only,
                                                            ext_hits=ext_hits)
                result_gene_lst.append(result_gene)
                if pos_result is not None:
                    result_pos_lst.append(pos_result)

            except:
                result_gene = pd.DataFrame({"Gene" : gene,
                                            "Uniprot_ID" : uniprot_id,
                                            "F" : np.nan,
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
        result_gene = add_nan_clust_cols(result_gene).drop(columns = ["Max_mut_pos", "Structure_max_pos"])
        result_gene = sort_cols(result_gene)
        if fragments == False:
            result_gene = result_gene.drop(columns=[col for col in ["F", "Mut_in_top_F", "Top_F"] if col in result_gene.columns])
        result_gene.to_csv(f"{output_dir}/{cohort}.3d_clustering_genes.csv", index=False)
        
    else:
        result_pos = pd.concat(result_pos_lst)
        result_pos["Cancer"] = cancer_type
        result_pos["Cohort"] = cohort
        result_pos.to_csv(f"{output_dir}/{cohort}.3d_clustering_pos.csv", index=False)
        # print("\n >> Wanted Pos>\n", result_pos.drop(columns=["Cancer", "Cohort"]))           #############################################################################################################################
        # Get gene global pval, qval, and clustering annotations
        #print(result_pos.head()) #######################à##############################################
        #print(result_pos.columns) ###################################################################################################################
        result_gene = get_final_gene_result(result_pos, result_gene, alpha)
        result_gene = result_gene.drop(columns = ["Max_mut_pos", "Structure_max_pos"]) 
        result_gene = sort_cols(result_gene) 
        if fragments == False:
            result_gene = result_gene.drop(columns=[col for col in ["F", "Mut_in_top_F", "Top_F"] if col in result_gene.columns])
        result_gene.to_csv(f"{output_dir}/{cohort}.3d_clustering_genes.csv", index=False)
        # print("\n >> Wanted Genes>\n", result_gene.drop(columns=["Uniprot_ID", "C_pos", "C_community", "Cancer", "Cohort"]))        #############################################################################################################################

if __name__ == "__main__":
    main()
 