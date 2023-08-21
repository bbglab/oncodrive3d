#!/usr/bin/env python

""" 
The module includes the the main function to run an HotMAPs-inspired 
method that uses the the canonical predicted structure stored in 
AlphaFold db to perform 3D-clustering of mutations using simulations 
and rank based comparison.

###################################### EXAMPLE USAGE ############################################

## CH

python3 main.py -i /workspace/projects/clustering_3d/evaluation/datasets/datasets_ch/input/maf/OTHER_WXS_CH_IMPACT_PANEL.in.maf \
-o /workspace/projects/clustering_3d/dev_testing/output/ \
-p /workspace/projects/clustering_3d/evaluation/datasets/datasets_ch/input/mut_profile/OTHER_WXS_CH_IMPACT_PANEL.mutrate.json \
-H 0 -t CH -C OTHER_WXS_CH_IMPACT_PANEL -e 1

## pCMAPs

python3 /workspace/projects/clustering_3d/clustering_3d/scripts/main.py  \
    -i /workspace/projects/clustering_3d/evaluation/datasets/input/maf/HARTWIG_WGS_PANCREAS_2020.in.maf  \
        -o /workspace/projects/clustering_3d/dev_testing/output \
            -p /workspace/projects/clustering_3d/evaluation/datasets/input/mut_profile/HARTWIG_WGS_PANCREAS_2020.mutrate.json \
                -s /workspace/projects/clustering_3d/clustering_3d/datasets_frag/seq_for_mut_prob.csv \
                    -c /workspace/projects/clustering_3d/clustering_3d/datasets_frag/prob_cmaps/ \
                        -d /workspace/projects/clustering_3d/clustering_3d/datasets_frag/confidence.csv \
                            -t BLCA -C HARTWIG_WGS_PANCREAS_2020 \
                                -e /workspace/projects/clustering_3d/clustering_3d/datasets_frag/pae \
                                    -u 30 -P 0.5

#################################################################################################
"""

import json
import logging
import os
import sys
from datetime import datetime

import click
import numpy as np
import pandas as pd
import daiquiri

from scripts import __logger_name__, __version__
from scripts.utils.clustering import clustering_3d_mp_wrapper
from scripts.utils.miss_mut_prob import (get_miss_mut_prob_dict,
                                         mut_rate_vec_to_dict)
from scripts.utils.pvalues import get_final_gene_result
from scripts.utils.utils import add_nan_clust_cols, parse_maf_input, sort_cols


DATE = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
FORMAT = "%(asctime)s - %(color)s%(levelname)s%(color_stop)s | %(name)s - %(color)s%(message)s%(color_stop)s"

def setup_logging(verbose: bool, output_path: str) -> None:
    """Set up logging facilities.

    :param verbose: verbosity (bool)
    :param output_path: str for log file
    """
    
    level = logging.DEBUG if verbose else logging.INFO

    formatter = daiquiri.formatter.ColorFormatter(fmt=FORMAT)

    daiquiri.setup(level=level, outputs=(
        daiquiri.output.Stream(formatter=formatter), 
        daiquiri.output.File(filename=output_path, formatter=formatter)
    ))

logger = daiquiri.getLogger(__logger_name__)


@click.group(context_settings={'help_option_names': ['-h', '--help']})
@click.version_option(__version__)
def oncodrive3D():
    """Oncrodrive3D: software for the identification of 3D-clustering of missense mutations for cancer driver genes detection."""
    pass

@oncodrive3D.command(context_settings=dict(help_option_names=['-h', '--help']),
               help="build datasets - required at the beginning") # CHANGE ACCORDINGLY
def build_datasets():
    """"function to build datasets"""

    pass


@oncodrive3D.command(context_settings=dict(help_option_names=['-h', '--help']),
               help="run analysis") # CHANGE ACCORDINGLY
@click.option("-i", "--input_maf_path", type=click.Path(exists=True), required=True, help="Path of the maf file used as input")
@click.option("-p", "--mut_profile_path", type=click.Path(exists=True), 
              help="Path to the mut profile (list of 96 floats) of the cohort (json)")
@click.option("-o", "--output_path", help="Path to output directory", type=str, required=True)
@click.option("-s", "--seq_df_path", type=click.Path(exists=True), 
              help="Path to the dataframe including DNA and protein seq of all gene/proteins (all AF predicted ones)")
@click.option("-c", "--cmap_path", type=click.Path(exists=True), 
              help="Path to the directory containing the contact map of each protein")
@click.option("-d", "--plddt_path", type=click.Path(exists=True), 
              help="Path to the pandas dataframe including the AF model confidence of all proteins")
@click.option("-e", "--pae_path", help="Path to the directory including the AF Predicted Aligned Error (PAE) .npy files")
@click.option("-n", "--n_iterations", help="Number of densities to be simulated", type=int, default=10000)
@click.option("-a", "--alpha", help="Significant threshold for the p-value of res and gene", type=float, default=0.01)
@click.option("-P", "--cmap_prob_thr", type=float, default=0.5,
              help="Threshold to define AAs contacts based on distance on predicted structure and PAE")
@click.option("-H", "--hits_only", help="If 1 returns only positions in clusters, if 0 returns all", is_flag=True)
@click.option("-f", "--no_fragments", help="Disable processing of fragmented (AF-F) proteins", is_flag=True)
@click.option("-u", "--num_cores", type=click.IntRange(min=1, max=os.cpu_count(), clamp=False), default=1,
              help="Set the number of cores to use in the computation")
@click.option("-S", "--seed", help="Set seed to ensure reproducible results", type=int)
@click.option("-v", "--verbose", help="Monitor the number of processed genes", is_flag=True)
@click.option("-t", "--cancer_type", help="Cancer type", type=str)
@click.option("-C", "--cohort", help="Name of the cohort", type=str)
def run(input_maf_path, 
         mut_profile_path, 
         output_path,
         seq_df_path,
         cmap_path,
         plddt_path,
         pae_path,
         n_iterations,
         alpha,
         cmap_prob_thr,
         hits_only,
         no_fragments,
         num_cores,
         seed,
         verbose,
         cancer_type,
         cohort):

    ## Initialize
    dir_path = os.path.abspath(os.path.dirname(__file__))
    if plddt_path is None:
        plddt_path = f"{dir_path}/../datasets/confidence.csv"
    if cmap_path is None:
        cmap_path = f"{dir_path}/../datasets/cmaps/"                                          
    if seq_df_path is None:
        seq_df_path = f"{dir_path}/../datasets/seq_for_mut_prob.csv"
    if cancer_type is None:
        cancer_type = np.nan
    if cohort is None:
        cohort = f"cohort_{DATE}"

    # Log
    if not os.path.exists(output_path):
        os.makedirs(os.path.join(output_path, 'logs/'))
        setup_logging(verbose,os.path.join(output_path,  'logs/', f'{cohort}-{DATE}.log'))
        logger.warning(f"Directory '{output_path}' does not exists")
        logger.debug(f"Creating output directory tree..")
    else: 
        setup_logging(verbose, os.path.join(output_path,  'logs/', f'{cohort}-{DATE}.log'))

    logger.info('Initializing Oncodrive3D...')
    logger.info(f'Oncodrive3D v{__version__}')

    logger.info(f"Starting 3D-clustering..")

    logger.info(f"Output directory: {output_path}")
    logger.info(f"Path to CMAPs: {cmap_path}")
    logger.info(f"Path to DNA sequences: {seq_df_path}")
    if pae_path is not None:
        logger.info(f"Path to PAE: {pae_path}")
    else:
        logger.info(f"Path to PAE: not defined, weighted average PAE of mutated volumes will not be calculated")
    logger.info(f"Path to pLDDT scores: {plddt_path}")

    path_prob = mut_profile_path if not None else "Not provided, uniform distribution will be used"

    logger.info(f"Path to mut profile: {path_prob}")
    logger.info(f"CPU cores: {num_cores}")
    logger.info(f"Iterations: {n_iterations}")
    logger.info(f"Significant level: {alpha}")
    logger.info(f"Probability threshold for CMAPs: {cmap_prob_thr}")
    logger.info(f"Output hits only: {bool(hits_only)}")
    logger.info(f"Disable fragments: {bool(no_fragments)}")
    logger.info(f"Cohort: {cohort}")
    logger.info(f"Cancer type: {cancer_type}")
    logger.info(f"Verbose: {bool(verbose)}")
    logger.info(f"Seed: {seed}")


    ## Load input and df of DNA sequences

    # MAF input
    data = parse_maf_input(input_maf_path, keep_samples_id=True)

    # Seq df for missense mut prob
    seq_df = pd.read_csv(seq_df_path)
    
    # Model confidences
    plddt_df = pd.read_csv(plddt_path, dtype={"Pos" : int, "Res" : str, "Confidence" : float, 
                                              "Uniprot_ID" : str, "AF_F" : str})


    ## Run

    result_np_gene_lst = []

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
    result_np_gene_lst.append(result_gene)   

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
    result_np_gene_lst.append(result_gene)
    
    # Filter on fragmented (AF-F) genes
    if no_fragments == True:
        # Return the fragmented genes as non processed output
        genes_frag = seq_df[seq_df.F.str.extract(r'(\d+)', expand=False).astype(int) > 1]
        genes_frag = genes_frag.Gene.reset_index(drop=True).values
        genes_frag_mut = genes_mut[[gene in genes_frag for gene in genes_mut.index]]
        genes_frag = genes_frag_mut.index.values
        genes_frag_id = [gene_to_uniprot_dict[gene] for gene in genes_frag]               
        result_gene = pd.DataFrame({"Gene" : genes_frag,
                                    "Uniprot_ID" : genes_frag_id,
                                    "F" : np.nan,
                                    "Mut_in_gene" : genes_frag_mut.values,
                                    "Max_mut_pos" : np.nan,
                                    "Structure_max_pos" : np.nan,
                                    "Status" : "Fragmented"})
        result_np_gene_lst.append(result_gene)
        # Filter out from genes to process and seq df
        genes_to_process = [gene for gene in genes_to_process if gene not in genes_frag]
        seq_df = seq_df[[gene in genes_to_process for gene in seq_df["Gene"]]].reset_index(drop=True)

    # Missense mut prob  
    if mut_profile_path is not None:
        # Compute dict from mut profile of the cohort and dna sequences
        mut_profile = json.load(open(mut_profile_path))
        logger.info(f"Computing missense mut probabilities..")
        if not isinstance(mut_profile, dict):
            mut_profile = mut_rate_vec_to_dict(mut_profile)
        miss_prob_dict = get_miss_mut_prob_dict(mut_rate_dict=mut_profile, seq_df=seq_df)
    else:
        miss_prob_dict = None

    # Run 3D-clustering
    if len(genes_to_process) > 0:
        logger.info(f"Performing 3D clustering on [{len(seq_df)}] proteins..")
        result_pos, result_gene = clustering_3d_mp_wrapper(genes_to_process, data, cmap_path, 
                                                           miss_prob_dict, gene_to_uniprot_dict, plddt_df,
                                                           num_cores, alpha=alpha, num_iteration=n_iterations, 
                                                           cmap_prob_thr=cmap_prob_thr, hits_only=hits_only, 
                                                           verbose=verbose, seed=seed, pae_path=pae_path)
        result_gene = pd.concat((result_gene, pd.concat(result_np_gene_lst)))
    else:
        result_gene = pd.concat(result_np_gene_lst)
        result_pos = None


    ## Save 
    logger.info(f"Saving to {output_path}")
    result_gene["Cancer"] = cancer_type
    result_gene["Cohort"] = cohort

    #if len(result_pos_lst) == 0:
    if result_pos is None:
        logger.warning(f"Did not processed any genes\n")
        result_gene = add_nan_clust_cols(result_gene).drop(columns = ["Max_mut_pos", "Structure_max_pos"])
        result_gene = sort_cols(result_gene)
        if no_fragments == True:
            result_gene = result_gene.drop(columns=[col for col in ["F", "Mut_in_top_F", "Top_F"] if col in result_gene.columns])
        result_gene.to_csv(f"{output_path}/{cohort}.3d_clustering_genes.csv", index=False)
        
    else:
        # Save res-level result
        result_pos["Cancer"] = cancer_type
        result_pos["Cohort"] = cohort
        result_pos.to_csv(f"{output_path}/{cohort}.3d_clustering_pos.csv", index=False)
   
        # Get gene global pval, qval, and clustering annotations and save gene-level result
        result_gene = get_final_gene_result(result_pos, result_gene, alpha)
        result_gene = result_gene.drop(columns = ["Max_mut_pos", "Structure_max_pos"]) 
        result_gene = sort_cols(result_gene) 
        if no_fragments == True:
            result_gene = result_gene.drop(columns=[col for col in ["F", "Mut_in_top_F", "Top_F"] if col in result_gene.columns])
        result_gene.to_csv(f"{output_path}/{cohort}.3d_clustering_genes.csv", index=False)

if __name__ == "__main__":
    oncodrive3D()
 