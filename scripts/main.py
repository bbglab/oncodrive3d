#!/usr/bin/env python

""" 
The module includes the the main function to run an HotMAPs-inspired 
method that uses the the canonical predicted structure stored in 
AlphaFold db to perform 3D-clustering of mutations using simulations 
and rank based comparison.

# =============
# EXAMPLE USAGE
# =============

# Build datasets

cd path/to/oncodrive3D
oncodrive3D build-datasets
oncodrive3D build-datasets -m -v -o datasets_mane

# Run 
            
oncodrive3D run \
    -i test/maf/TCGA_WXS_ACC.in.maf  \
        -p test/mut_profile/TCGA_WXS_ACC.mutrate.json \
            -o test/results
                  

singularity exec /workspace/projects/clustering_3d/clustering_3d/build/containers/oncodrive3d_231205.sif oncodrive3D run -i /workspace/projects/clustering_3d/o3d_analysys/datasets/input/normal/kidney_pilot/all_mutations.all_samples.tsv -d /workspace/projects/clustering_3d/clustering_3d/datasets_normal -m /workspace/projects/clustering_3d/o3d_analysys/datasets/input/normal/kidney_pilot/mutability_kidney.json -o <your_output_dir> -C kidney_normal -v -s 128 -t kidney 
singularity exec /workspace/projects/clustering_3d/clustering_3d/build/containers/oncodrive3d_231205.sif oncodrive3D run -i /workspace/projects/clustering_3d/o3d_analysys/datasets/input/normal/kidney_pilot/all_mutations.all_samples.tsv -d /workspace/projects/clustering_3d/clustering_3d/datasets_normal -m /workspace/projects/clustering_3d/o3d_analysys/datasets/input/normal/kidney_pilot/mutability_kidney.json -o /workspace/projects/clustering_3d/o3d_analysys/datasets/output/normal/o3d_output/kidney_231205_test -C kidney_normal -v -s 128 -t kidney 
singularity exec oncodrive3D plot -i /workspace/projects/clustering_3d/o3d_analysys/datasets/output/normal/o3d_output/kidney_231205_test -c kidney_normal

oncodrive3D run -i /workspace/projects/clustering_3d/o3d_analysys/datasets/input/normal/kidney_pilot/all_mutations.all_samples.tsv -d /workspace/nobackup/scratch/oncodrive3d/datasets -m /workspace/projects/clustering_3d/o3d_analysys/datasets/input/normal/kidney_pilot/mutability_kidney.json -o /workspace/projects/clustering_3d/o3d_analysys/datasets/output/normal/o3d_output/kidney_240404 -C kidney_normal -v -s 128 -t kidney 


-- Old & nEw comparison

- Old
oncodrive3D run -i /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer/maf/PCAWG_WGS_ESO_ADENOCA.in.maf -p /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer/mut_profile/PCAWG_WGS_ESO_ADENOCA.mutrate.json -d /workspace/nobackup/scratch/oncodrive3d/datasets_last_real -C PCAWG_WGS_ESO_ADENOCA -o /workspace/projects/clustering_3d/dev_testing/result/o3d/PCAWG_WGS_ESO_ADENOCA -s 128 -c 10

- New
oncodrive3D run -i /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer_202404/maf/PCAWG_WGS_ESO_ADENOCA.in.maf -p /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer_202404/mut_profile/PCAWG_WGS_ESO_ADENOCA.sig.json -d /workspace/nobackup/scratch/oncodrive3d/datasets_last_real -C PCAWG_WGS_ESO_ADENOCA -o /workspace/projects/clustering_3d/dev_testing/result/o3d/PCAWG_WGS_ESO_ADENOCA_new -s 128 -c 10

- New vep output as input
oncodrive3D run -i /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer_202404/vep/PCAWG_WGS_ESO_ADENOCA.vep.tsv.gz -p /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer_202404/mut_profile/PCAWG_WGS_ESO_ADENOCA.sig.json -d /workspace/nobackup/scratch/oncodrive3d/datasets_last_real -C PCAWG_WGS_ESO_ADENOCA -o /workspace/projects/clustering_3d/dev_testing/result/o3d/PCAWG_WGS_ESO_ADENOCA_new_vep -s 128 -c 10 --o3d_transcripts
oncodrive3D run -i /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer_202404/vep/TCGA_WXS_BLCA.vep.tsv.gz -p /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer_202404/mut_profile/TCGA_WXS_BLCA.sig.json -d /workspace/nobackup/scratch/oncodrive3d/datasets_last_real -C TCGA_WXS_BLCA -o /workspace/projects/clustering_3d/dev_testing/result/o3d/TCGA_WXS_BLCA_new -s 128 -c 10 --o3d_transcripts --use_input_symbols -v
oncodrive3D run -i /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer_202404/vep/CBIOP_WXS_ACY_2019.vep.tsv.gz -p /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer_202404/mut_profile/CBIOP_WXS_ACY_2019.sig.json -d /workspace/nobackup/scratch/oncodrive3d/datasets_240506 -C CBIOP_WXS_ACY_2019 -o /workspace/projects/clustering_3d/dev_testing/result/o3d/test_240612 -s 26 -c 10 -v --o3d_transcripts --use_input_symbols -v

# New vep output as input MANE
oncodrive3D run -i /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer_202404/vep/PCAWG_WGS_ESO_ADENOCA.vep.tsv.gz -p /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer_202404/mut_profile/PCAWG_WGS_ESO_ADENOCA.sig.json -d /workspace/nobackup/scratch/oncodrive3d/datasets_mane_last_real -C PCAWG_WGS_ESO_ADENOCA -o /workspace/projects/clustering_3d/dev_testing/result/o3d/PCAWG_WGS_ESO_ADENOCA_new_mane_vep -s 128 -c 10 --o3d_transcripts --use_input_symbols -v --mane
oncodrive3D run -i /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer_202404/vep/TCGA_WXS_BLCA.vep.tsv.gz -p /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer_202404/mut_profile/TCGA_WXS_BLCA.sig.json -d /workspace/nobackup/scratch/oncodrive3d/datasets_mane_last_real -C TCGA_WXS_BLCA -o /workspace/projects/clustering_3d/dev_testing/result/o3d/TCGA_WXS_BLCA_new_mane -s 128 -c 10 --o3d_transcripts --use_input_symbols -v --mane


# Plot

# oncodrive3D plot --output_tsv --non_significant -r kidney_normal -g /workspace/projects/clustering_3d/o3d_analysys/datasets/output/cancer/o3d_output/run_ref_trinucl/results/TCGA_WXS_BLCA.3d_clustering_genes.tsv -p /workspace/projects/clustering_3d/o3d_analysys/datasets/output/cancer/o3d_output/run_ref_trinucl/results/TCGA_WXS_BLCA.3d_clustering_pos.tsv -i /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer/maf/TCGA_WXS_BLCA.in.maf -o /workspace/projects/clustering_3d/o3d_analysys/datasets/output/cancer/o3d_output/run_ref_trinucl/plots -m /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer/mut_profile/TCGA_WXS_BLCA.mutrate.json -d /workspace/projects/clustering_3d/clustering_3d/datasets -a /workspace/projects/clustering_3d/o3d_analysys/datasets/annotations
"""


# =============================================================================
# TODO: allow procesing without tumor sample info
# TODO: handle inf of the score (e.g., decimal python package)
# TODO: fix bug in requirement.txt (bgreference must be installed after setup.py)
# TODO: add script to generate conf and mutability file
# TODO: test run on normal tissue
# TODO: add filter (and logs) for mutated genes without exons coordinate when mutability is provided

# TODO: change progressbar to tqdm in run scripts
# TODO: change output names?
# TODO: change repo name
# TODO: fix doc
# TODO: update doc for normal tissue application
# TODO: suppress verbosity of multi-threading download of structures
# =============================================================================


import json
import os

import click
import daiquiri
import numpy as np
import pandas as pd

from scripts import __logger_name__, __version__
from scripts.datasets.build_datasets import build
from scripts.globals import DATE, setup_logging_decorator, startup_message
from scripts.run.clustering import clustering_3d_mp_wrapper
from scripts.run.miss_mut_prob import get_miss_mut_prob_dict, mut_rate_vec_to_dict
from scripts.run.pvalues import get_final_gene_result
from scripts.run.utils import get_gene_entry, add_nan_clust_cols, parse_maf_input, sort_cols, empty_result_pos
from scripts.plotting.build_annotations import get_annotations
from scripts.plotting.plot import generate_plots, generate_comparative_plots
from scripts.plotting.utils import init_plot_pars, init_comp_plot_pars, parse_lst_tracks
from scripts.run.mutability import init_mutabilities_module

logger = daiquiri.getLogger(__logger_name__)



@click.group(context_settings={'help_option_names': ['-h', '--help']})
@click.version_option(__version__)
def oncodrive3D():
    """
    Oncodrive3D: software for the identification of 3D-clustering 
    of missense mutations for cancer driver genes detection.
    """
    pass


# =============================================================================
#                               BUILD DATASETS
# =============================================================================

@oncodrive3D.command(context_settings=dict(help_option_names=['-h', '--help']),
                     help="Build datasets - Required once after installation.")
@click.option("-o", "--output_dir", 
              help="Directory where to save the files", type=str, default='datasets')
@click.option("-s", "--organism", type=click.Choice(["Homo sapiens", 'human', "Mus musculus", 'mouse']), 
              help="Organism name", default="Homo sapiens")
@click.option("-m", "--mane", 
              help="Use structures predicted from MANE Select transcripts (Homo sapiens only)", is_flag=True)
@click.option("-M", "--mane_version", default=1.3, 
              help="Version of the MANE Select release from NCBI")
@click.option("-d", "--distance_threshold", type=click.INT, default=10,
              help="Distance threshold (Å) to define contact between amino acids")
@click.option("-c", "--cores", type=click.IntRange(min=1, max=len(os.sched_getaffinity(0)), clamp=False), default=len(os.sched_getaffinity(0)),
              help="Number of cores to use in the computation")
@click.option("--af_version", type=click.IntRange(min=1, max=4, clamp=False), default=4,
              help="Version of AlphaFold 2 predictions")
@click.option("-y", "--yes", 
              help="No interaction", is_flag=True)
@click.option("-v", "--verbose", 
              help="Verbose", is_flag=True)
@setup_logging_decorator
def build_datasets(output_dir,
                   organism,
                   mane,
                   distance_threshold,
                   cores, 
                   af_version,
                   mane_version,
                   yes,
                   verbose):
    """"Build datasets necessary to run Oncodrive3D."""
    startup_message(__version__, "Initializing building datasets..")

    logger.info(f"Current working directory: {os.getcwd()}")
    logger.info(f"Build folder path: {output_dir}")
    logger.info(f"Organism: {organism}")
    logger.info(f"MANE Select: {mane}")
    logger.info(f"Distance threshold: {distance_threshold}Å")
    logger.info(f"CPU cores: {cores}")
    logger.info(f"AlphaFold version: {af_version}")
    logger.info(f"MANE version: {mane_version}")
    logger.info(f"Verbose: {verbose}")
    logger.info(f'Log path: {os.path.join(output_dir, "log")}')
    logger.info("")

    build(output_dir,
          organism,
          mane,
          distance_threshold,
          cores,
          af_version,
          mane_version)



# =============================================================================
#                                     RUN
# =============================================================================

@oncodrive3D.command(context_settings=dict(help_option_names=['-h', '--help']),
                     help="Run 3D-clustering analysis.")
@click.option("-i", "--input_maf_path", type=click.Path(exists=True), required=True,
              help="Path of the MAF file (or direct VEP output) used as input")
@click.option("-p", "--mut_profile_path", type=click.Path(exists=True), 
              help="Path of the mutation profile (192 trinucleotide contexts) used as optional input")
@click.option("-m", "--mutability_config_path", type=click.Path(exists=True), 
              help="Path of the config file with information on mutability")
@click.option("-o", "--output_dir", type=str, default='results', 
              help="Path to output directory")
@click.option("-d", "--data_dir", type=click.Path(exists=True), default = os.path.join('datasets'), 
              help="Path to datasets")
@click.option("-n", "--n_iterations", type=int, default=10000, 
              help="Number of densities to be simulated")
@click.option("-a", "--alpha", type=float, default=0.01, 
              help="Significant threshold for the p-value of res and gene")
@click.option("-P", "--cmap_prob_thr", type=float, default=0.5,
              help="Threshold to define AAs contacts based on distance on predicted structure and PAE")
@click.option("-c", "--cores", type=click.IntRange(min=1, max=len(os.sched_getaffinity(0)), clamp=False), default=len(os.sched_getaffinity(0)),
              help="Set the number of cores to use in the computation")
@click.option("-s", "--seed", type=int,
              help="Set seed to ensure reproducible results")
@click.option("-v", "--verbose", 
              help="Verbose", is_flag=True)
@click.option("-t", "--cancer_type", 
              help="Cancer type", type=str)
@click.option("-C", "--cohort", 
              help="Name of the cohort", type=str)
@click.option("--no_fragments", is_flag=True, 
              help="Disable processing of fragmented (AF-F) proteins")
@click.option("--only_processed", is_flag=True,
              help="Include only processed genes in the output")
@click.option("--thr_mapping_issue", type=float, default=0.1,
              help="Threshold to filter out genes by the ratio of mutations with mapping issue (out of structure, WT AA mismatch, zero prob to mutate)")
@click.option("--o3d_transcripts", is_flag=True,
              help="Filter mutations by keeping transcripts included in Oncodrive3D built sequence dataframe. Only if input file (--i) is a raw VEP output")
@click.option("--use_input_symbols", is_flag=True,
              help="Update HUGO symbols in Oncodrive3D built datasets by using input file entries. Only if input file (--i) is a raw VEP output")
@click.option("--mane", is_flag=True,
              help="If multiple structures are associated to the same HUGO symbol in the input file, use the MANE ones.")
@setup_logging_decorator
def run(input_maf_path,
        mut_profile_path,
        mutability_config_path,
        output_dir,
        data_dir,
        n_iterations,
        alpha,
        cmap_prob_thr,
        cores,
        seed,
        verbose,
        cancer_type,
        cohort,
        no_fragments,
        only_processed,
        thr_mapping_issue,
        o3d_transcripts,
        use_input_symbols,
        mane):
    """Run Oncodrive3D."""

    ## Initialize

    plddt_path = os.path.join(data_dir, "confidence.tsv")
    cmap_path = os.path.join(data_dir, "prob_cmaps")
    seq_df_path = os.path.join(data_dir, "seq_for_mut_prob.tsv")
    pae_path = os.path.join(data_dir, "pae")
    cancer_type = cancer_type if cancer_type else np.nan
    cohort = cohort if cohort else f"cohort_{DATE}"
    path_prob = mut_profile_path if mut_profile_path else "Not provided, mutabilities will be used" if mutability_config_path else "Not provided, uniform distribution will be used"
    path_mutability_config = mutability_config_path if mutability_config_path else "Not provided, mutabilities will not be used"

    # Log
    startup_message(__version__, "Initializing analysis..")

    logger.info(f"Input MAF: {input_maf_path}")
    logger.info(f"Input mut profile: {path_prob}")
    logger.info(f"Input mutability config: {path_mutability_config}")
    logger.info(f"Build directory: {data_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.debug(f"Path to CMAPs: {cmap_path}")
    logger.debug(f"Path to DNA sequences: {seq_df_path}")
    logger.debug(f"Path to PAE: {pae_path}")
    logger.debug(f"Path to pLDDT scores: {plddt_path}")
    logger.info(f"CPU cores: {cores}")
    logger.info(f"Iterations: {n_iterations}")
    logger.info(f"Significant level: {alpha}")
    logger.info(f"Probability threshold for CMAPs: {cmap_prob_thr}")
    logger.info(f"Cohort: {cohort}")
    logger.info(f"Cancer type: {cancer_type}")
    logger.info(f"Disable fragments: {bool(no_fragments)}")
    logger.info(f"Output only processed genes: {bool(only_processed)}")
    logger.info(f"Ratio threshold mutations with mapping issue: {thr_mapping_issue}")
    logger.info(f"Seed: {seed}")
    logger.info(f"Filter input by Oncodrive3D transcripts (only if VEP output is used as input): {o3d_transcripts}")
    logger.info(f"Use HUGO symbols of input file (only if VEP output is used as input): {use_input_symbols}")
    logger.info(f"Prioritize MANE transcripts when using input HUGO symbols: {mane}")
    logger.info(f"Verbose: {bool(verbose)}")
    logger.info(f'Log path: {os.path.join(output_dir, "log")}')
    logger.info("")


    # Load
    #=====

    seq_df = pd.read_csv(seq_df_path, sep="\t")
    data, seq_df = parse_maf_input(input_maf_path, 
                                   seq_df, 
                                   use_o3d_transcripts=o3d_transcripts,
                                   use_input_symbols=use_input_symbols, 
                                   mane=mane)
    
    if len(data) > 0:


        # Run
        #====

        # Get genes with enough mut
        result_np_gene_lst = []
        genes = data.groupby("Gene").apply(len)
        genes_mut = genes[genes >= 2]
        genes_no_mut = genes[genes < 2].index

        if len(genes_no_mut) > 0:
            logger.debug(f"Detected [{len(genes_no_mut)}] genes without enough mutations: Skipping..")
            result_gene = pd.DataFrame({"Gene" : genes_no_mut,
                                        "Uniprot_ID" : np.nan,
                                        "F" : np.nan,
                                        "Mut_in_gene" : 1,
                                        "Ratio_not_in_structure" : np.nan,
                                        "Ratio_WT_mismatch" : np.nan,
                                        "Mut_zero_mut_prob" : np.nan,
                                        "Pos_zero_mut_prob" : np.nan,
                                        "Transcript_ID" : get_gene_entry(data, genes_no_mut, "Transcript_ID"),
                                        "O3D_transcript_ID" : get_gene_entry(data, genes_no_mut, "O3D_transcript_ID"),
                                        "Transcript_status" : get_gene_entry(data, genes_no_mut, "Transcript_status"),
                                        "Status" : "No_mut"})
            result_np_gene_lst.append(result_gene)

        # Seq df for metadata info
        metadata_cols = [col for col in ["Gene", "HGNC_ID", "Ens_Gene_ID", "Ens_Transcr_ID", "Refseq_prot", "Uniprot_ID", "F"] if col in seq_df.columns]
        metadata_mapping_cols = [col for col in ["Seq", "Chr", "Reverse_strand", "Exons_coord", "Seq_dna", "Tri_context", "Reference_info"] if col in seq_df.columns]
        seq_df_all = seq_df[seq_df["Gene"].isin(genes.index)].copy()

        # Get genes with corresponding Uniprot-ID mapping        
        gene_to_uniprot_dict = {gene : uni_id for gene, uni_id in seq_df[["Gene", "Uniprot_ID"]].drop_duplicates().values}
        genes_to_process = [gene for gene in genes_mut.index if gene in gene_to_uniprot_dict.keys()]
        seq_df = seq_df[seq_df["Gene"].isin(genes_to_process)].reset_index(drop=True)
        genes_no_mapping = genes[[gene in genes_mut.index and gene not in gene_to_uniprot_dict.keys() for gene in genes.index]]
        if len(genes_no_mapping) > 0:
            logger.debug(f"Detected [{len(genes_no_mapping)}] genes without IDs mapping: Skipping..")
            result_gene = pd.DataFrame({"Gene" : genes_no_mapping.index,
                                        "Uniprot_ID" : np.nan,
                                        "F" : np.nan,
                                        "Mut_in_gene" : genes_no_mapping.values,
                                        "Ratio_not_in_structure" : np.nan,
                                        "Ratio_WT_mismatch" : np.nan,
                                        "Mut_zero_mut_prob" : np.nan,
                                        "Pos_zero_mut_prob" : np.nan,
                                        "Transcript_ID" : get_gene_entry(data, genes_no_mapping.index, "Transcript_ID"),
                                        "O3D_transcript_ID" : get_gene_entry(data, genes_no_mapping.index, "O3D_transcript_ID"),
                                        "Transcript_status" : get_gene_entry(data, genes_no_mapping.index, "Transcript_status"),
                                        "Status" : "No_ID_mapping"})
            result_np_gene_lst.append(result_gene)
        
        # Filter on fragmented (AF-F) genes
        if no_fragments:
            # Return the fragmented genes as non processed output
            genes_frag = seq_df[seq_df.F.str.extract(r'(\d+)', expand=False).astype(int) > 1]
            genes_frag = genes_frag.Gene.reset_index(drop=True).values
            genes_frag_mut = genes_mut[[gene in genes_frag for gene in genes_mut.index]]
            genes_frag = genes_frag_mut.index.values
            if len(genes_frag) > 0:
                logger.debug(f"Detected [{len(genes_frag)}] fragmented genes with disabled fragments processing: Skipping..")
                result_gene = pd.DataFrame({"Gene" : genes_frag,
                                            "Uniprot_ID" : np.nan, 
                                            "F" : np.nan,
                                            "Mut_in_gene" : genes_frag_mut.values,
                                            "Ratio_not_in_structure" : np.nan,
                                            "Ratio_WT_mismatch" : np.nan,
                                            "Mut_zero_mut_prob" : np.nan,
                                            "Pos_zero_mut_prob" : np.nan,
                                            "Transcript_ID" : get_gene_entry(data, genes_frag, "Transcript_ID"),
                                            "O3D_transcript_ID" : get_gene_entry(data, genes_frag, "O3D_transcript_ID"),
                                            "Transcript_status" : get_gene_entry(data, genes_frag, "Transcript_status"),
                                            "Status" : "Fragmented"})
                result_np_gene_lst.append(result_gene)
                # Filter out from genes to process and seq df
                genes_to_process = [gene for gene in genes_to_process if gene not in genes_frag]
                seq_df = seq_df[seq_df["Gene"].isin(genes_to_process)].reset_index(drop=True)
                
        # Filter on start-loss mutations
        start_mut_ix = data["Pos"] == 1
        start_mut = sum(start_mut_ix)
        if start_mut > 0:
            genes_start_mut = list(data[start_mut_ix].Gene.unique())
            data = data[~start_mut_ix]
            logger.warning(f"Detected {start_mut} start-loss mutations in {len(genes_start_mut)} genes {genes_start_mut}: Filtering mutations..")


        # Missense mut prob
        #==================
        
        # Using mutabilities if provided
        if mutability_config_path is not None:
            logger.info("Computing missense mut probabilities using mutabilities..")
            mutab_config = json.load(open(mutability_config_path, encoding="utf-8"))
            init_mutabilities_module(mutab_config)
            seq_df = seq_df[seq_df["Reference_info"] == 1]   
            seq_df['Exons_coord'] = seq_df['Exons_coord'].apply(eval)  
            genes_to_process = [gene for gene in genes_to_process if gene in seq_df["Gene"].unique()]
            genes_not_mutability = [gene for gene in genes_to_process if gene not in seq_df["Gene"].unique()]
            miss_prob_dict = get_miss_mut_prob_dict(mut_rate_dict=None, seq_df=seq_df,
                                                    mutability=True, mutability_config=mutab_config)
            
            if len(genes_not_mutability) > 0:   
                logger.debug(f"Detected [{len(genes_not_mutability)}] genes without mutability information: Skipping..")
                result_gene = pd.DataFrame({"Gene" : genes_not_mutability,
                                            "Uniprot_ID" : np.nan,
                                            "F" : np.nan,
                                            "Mut_in_gene" : np.nan,
                                            "Ratio_not_in_structure" : np.nan,
                                            "Ratio_WT_mismatch" : np.nan,
                                            "Mut_zero_mut_prob" : np.nan,
                                            "Pos_zero_mut_prob" : np.nan,
                                            "Transcript_ID" : get_gene_entry(data, genes_not_mutability, "Transcript_ID"),
                                            "O3D_transcript_ID" : get_gene_entry(data, genes_not_mutability, "O3D_transcript_ID"),
                                            "Transcript_status" : get_gene_entry(data, genes_not_mutability, "Transcript_status"),
                                            "Status" : "No_mutability"})
                result_np_gene_lst.append(result_gene)
                
        # Using mutational profiles
        elif mut_profile_path is not None:
            # Compute dict from mut profile of the cohort and dna sequences
            mut_profile = json.load(open(mut_profile_path, encoding="utf-8"))
            logger.info("Computing missense mut probabilities..")
            if not isinstance(mut_profile, dict):
                mut_profile = mut_rate_vec_to_dict(mut_profile)
            miss_prob_dict = get_miss_mut_prob_dict(mut_rate_dict=mut_profile, seq_df=seq_df)
        else:
            logger.warning("Mutation profile not provided: Uniform distribution will be used for scoring and simulations.")
            miss_prob_dict = None


        # Run 3D-clustering
        #==================
        
        if len(result_np_gene_lst):
            result_np_gene = pd.concat(result_np_gene_lst)
            result_np_gene["Uniprot_ID"] = [gene_to_uniprot_dict[gene] if gene in gene_to_uniprot_dict.keys() else np.nan for gene in result_np_gene["Gene"].values]
        if len(genes_to_process) > 0:
            logger.info(f"Performing 3D-clustering on [{len(seq_df)}] proteins..")
            seq_df = seq_df[["Gene", "Uniprot_ID", "F", "Seq"]]
            plddt_df = pd.read_csv(plddt_path, sep="\t", usecols=["Pos", "Confidence", "Uniprot_ID"], dtype={"Pos" : np.int32,
                                                                                                             "Confidence" : np.float32, 
                                                                                                             "Uniprot_ID" : "object"})  
            
            result_pos, result_gene = clustering_3d_mp_wrapper(genes=genes_to_process,
                                                               data=data,
                                                               cmap_path=cmap_path,
                                                               miss_prob_dict=miss_prob_dict,
                                                               seq_df=seq_df,
                                                               plddt_df=plddt_df,
                                                               num_cores=cores,
                                                               alpha=alpha,
                                                               num_iteration=n_iterations,
                                                               cmap_prob_thr=cmap_prob_thr,
                                                               seed=seed,
                                                               pae_path=pae_path,
                                                               thr_mapping_issue=thr_mapping_issue)
            if result_np_gene_lst:
                result_gene = pd.concat((result_gene, result_np_gene))
        else:
            result_gene = result_np_gene
            result_pos = None


        # Save
        #=====
        
        if not os.path.exists(output_dir):
            os.makedirs(os.path.join(output_dir))
            logger.warning(f"Directory '{output_dir}' does not exists: Creating..")

        result_gene["Cancer"] = cancer_type
        result_gene["Cohort"] = cohort
        output_path_pos = os.path.join(output_dir, f"{cohort}.3d_clustering_pos.csv")
        output_path_genes = os.path.join(output_dir, f"{cohort}.3d_clustering_genes.csv")
        
        # Save processed seq_df and input files
        seq_df_output = os.path.join(output_dir, f"{cohort}.seq_df.processed.tsv")
        input_mut_output = os.path.join(output_dir, f"{cohort}.mutations.processed.tsv")
        input_prob_output = os.path.join(output_dir, f"{cohort}.miss_prob.processed.json")
        logger.info(f"Saving {seq_df_output}")
        seq_df_all[metadata_cols + metadata_mapping_cols].to_csv(seq_df_output, sep="\t", index=False)
        logger.info(f"Saving {input_mut_output}")
        data.to_csv(input_mut_output, sep="\t", index=False)
        logger.info(f"Saving {input_prob_output}")
        with open(input_prob_output, "w") as json_file:
            json.dump(miss_prob_dict, json_file)
        
        # Add extra metadata
        result_gene = result_gene.merge(seq_df_all[metadata_cols], on=["Gene", "Uniprot_ID"], how="left")

        if only_processed:
            result_gene = result_gene[result_gene["Status"] == "Processed"]

        if result_pos is None:
            # Save gene-level result and empty res-level result
            logger.warning("Did not processed any genes!")
            result_gene = add_nan_clust_cols(result_gene)
            result_gene = sort_cols(result_gene)
            if no_fragments:
                result_gene = result_gene.drop(columns=[col for col in ["F", "Mut_in_top_F", "Top_F"] if col in result_gene.columns])
            empty_result_pos().to_csv(output_path_pos, index=False)
            result_gene.to_csv(output_path_genes, index=False)

            logger.info(f"Saving (empty) {output_path_pos}")
            logger.info(f"Saving {output_path_genes}")

        else:
            # Save res-level result
            result_pos["Cancer"] = cancer_type
            result_pos["Cohort"] = cohort
            result_pos.to_csv(output_path_pos, index=False)

            # Get gene global pval, qval, and clustering annotations and save gene-level result
            result_gene = get_final_gene_result(result_pos, result_gene, alpha)
            result_gene = result_gene
            result_gene = sort_cols(result_gene) 
            if no_fragments:
                result_gene = result_gene.drop(columns=[col for col in ["F", "Mut_in_top_F", "Top_F"] if col in result_gene.columns])
            with np.printoptions(linewidth=10000):
                result_gene.to_csv(output_path_genes, index=False)

            logger.info(f"Saving {output_path_pos}")
            logger.info(f"Saving {output_path_genes}")

        logger.info("3D-clustering analysis completed!")

    else:
        logger.warning("No missense mutations were found in the input MAF. Consider checking your data: the field 'Variant_Classification' should include either 'Missense_Mutation' or 'missense_variant'")
        
               

# =============================================================================
#                              BUILD ANNOTATIONS
# =============================================================================

# Example:
# oncodrive3D build-annotations -o annotations_final -v -d /workspace/projects/clustering_3d/clustering_3d/datasets

# TODO: maybe use as input the path to datasets, then retrieve the structure from there.

@oncodrive3D.command(context_settings=dict(help_option_names=['-h', '--help']),
               help="Build annotations - Required (once) only to plot annotations.")
@click.option("-d", "--data_dir", help="Path to datasets", type=str, required=True)
@click.option("-o", "--output_dir", help="Path to dir where to store annotations", type=str, default="annotations")
@click.option("-g", "--ddg_dir", help="Path to custom ddG predictions", type=str)
#@click.option("-S", "--path_pdb_tool_sif", help="Path to the PDB_Tool SIF", type=str, required=True) 
@click.option("-s", "--organism", type=click.Choice(["Homo sapiens", 'human', "Mus musculus", 'mouse']), help="Organism name", default="Homo sapiens")
@click.option("-c", "--cores", type=click.IntRange(min=1, max=len(os.sched_getaffinity(0)), clamp=False), default=len(os.sched_getaffinity(0)),
              help="Number of cores to use in the computation")
@click.option("-y", "--yes", help="No interaction", is_flag=True)
@click.option("-v", "--verbose", help="Verbose", is_flag=True)
@setup_logging_decorator
def build_annotations(data_dir,
                      output_dir,
                      ddg_dir,
                      #path_pdb_tool_sif,
                      organism,
                      cores,
                      yes,
                      verbose):
    """
    Build datasets to plot protein annotations.
    
    Example: oncodrive3D build-annotations -o /workspace/projects/clustering_3d/o3d_analysys/datasets/annotations -v
    """

    startup_message(__version__, "Initializing building annotations..")
    
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Path to datasets: {data_dir}")
    logger.info(f"Path to custom ddG predictions: {ddg_dir}")
    #logger.info(f"Path to PDB_Tool SIF: {path_pdb_tool_sif}")
    logger.info(f"Organism: {organism}")
    logger.info(f"Cores: {cores}")
    logger.info(f"Verbose: {bool(verbose)}")
    logger.info(f'Log path: {os.path.join(output_dir, "log")}')
    logger.info("")

    get_annotations(data_dir, 
                    output_dir, 
                    ddg_dir,
                    #path_pdb_tool_sif,
                    organism,
                    cores)



# =============================================================================
#                                    PLOT
# =============================================================================

# Example:
# oncodrive3D plot -i /workspace/projects/clustering_3d/o3d_analysys/datasets/output/normal/o3d_output/2024/ -d /workspace/nobackup/scratch/oncodrive3d/datasets_240506/ -a /workspace/nobackup/scratch/oncodrive3d/annotations_240506/ -o /workspace/projects/clustering_3d/dev_testing/result/plots -c all_samples_raw --title "Bladder normal tissue" --non_significant --output_tsv -v --summary_alpha 0.3

# TO DO
# - Add a check for lst hratios that match with lst tracks and annotations
# - Add a check for the names of tracks and annotations

@oncodrive3D.command(context_settings=dict(help_option_names=['-h', '--help']),
               help="Generate plots for a quick interpretation of the 3D-clustering analysis.")
@click.option("-g", "--gene_result_path", type=click.Path(exists=True), required=True,
              help="Path to genes-level O3D result")
@click.option("-p", "--pos_result_path", type=click.Path(exists=True), required=True,
              help="Path to positions-level O3D result")
@click.option("-i", "--maf_path", type=click.Path(exists=True), required=True,
              help="Path to input mutations file")
@click.option("-m", "--miss_prob_path", type=click.Path(exists=True), required=True,
              help="Path to missense mutations probability dictionary")
@click.option("-s", "--seq_df_path", type=click.Path(exists=True), required=True,
              help="Path to dataframe of sequences")
@click.option("-d", "--datasets_dir", type=click.Path(exists=True), required=True,
              help="Path to datasets directory")
@click.option("-a", "--annotations_dir", type=click.Path(exists=True), required=True, 
              help="Path to annotations directory")
@click.option("-o", "--output_dir", default="./",
              help="Path to output directory where to save plots")
@click.option("-c", "--cohort", 
              help="Cohort name", type=str, required=True)
@click.option("--title",                                                                              ## Might be redundant
              help="Plot title", type=str)
@click.option("--maf_for_nonmiss_path", type=click.Path(exists=True), 
              help="Path to input mutations file including non-missense mutations")
@click.option("--lst_summary_tracks", type=str,
              help="List of tracks to be included in the summary plot (e.g., score,miss_count,clusters)", 
              default="score,miss_count,res_count,res_ratio,clusters")
@click.option("--lst_summary_hratios", type=str,
              help="List of float to define horizontal ratio of each track of the summary plot") 
@click.option("--lst_gene_tracks", type=str,
              help="List of tracks to be included in the gene plots (e.g., miss_count,miss_prob,score)",
              default="miss_count,miss_prob,score,clusters,ddg,disorder,pacc,ptm,site,sse,pfam,prosite,membrane,motif")
@click.option("--lst_gene_hratios", type=str,
              help="List of floats to define horizontal ratio of each track of the gene plot") 
@click.option("--summary_fsize_x", help="Figure size x-axis for summary plots (to be dynamically increased)", type=float, default=0.5)
@click.option("--summary_fsize_y", help="Figure size y-axis for summary plots", type=int, default=8)
@click.option("--gene_fsize_x", help="Figure size x-axis for gene plots", type=int, default=24)
@click.option("--gene_fsize_y", help="Figure size y-axis for gene plots", type=int, default=12)
@click.option("--summary_alpha", help="Alpha value for score track in summary plot", type=float, default=0.7)
@click.option("--dist_thr", help="Threshold of ratios to avoid clashing feature names (e.g., domains and motifs)", type=float, default=0.1)
@click.option("--genes", help="List of genes to be analysed in the report (e.g., --genes TP53,KRAS,PIK3CA)", type=str)
@click.option("--n_genes", help="Top number of genes to be included in the plots", type=int, default=30)
@click.option("--volcano_top_n", help="Top associations to annotated in volcano plot", type=int, default=15)
@click.option("--volcano_fsize_x", help="Figure size x-axis for volcano plot", type=int, default=10)
@click.option("--volcano_fsize_y", help="Figure size y-axis for volcano plot", type=int, default=6)
@click.option("--volcano_subplots_fsize_x", help="Figure size x-axis for volcano subplots", type=int, default=15)
@click.option("--volcano_subplots_fsize_y", help="Figure size y-axis for volcano subplots", type=int, default=10)
@click.option("--log_odds_fsize_x", help="Figure size x-axis for log odds plot", type=int, default=20)
@click.option("--log_odds_fsize_y", help="Figure size y-axis for log odds plot", type=int, default=5.5)
@click.option("--output_csv", help="Output csv file including annotated Oncodrive3D result", is_flag=True)
@click.option("--output_all_pos", help="Include all position (including non-mutated ones) in the Oncodrive3D enriched result", is_flag=True)
@click.option("-v", "--verbose", help="Verbose", is_flag=True)
@setup_logging_decorator
def plot(gene_result_path,
         pos_result_path,
         maf_path,
         miss_prob_path,
         seq_df_path,
         datasets_dir,
         annotations_dir,
         output_dir,
         cohort,
         title,
         maf_for_nonmiss_path,
         lst_summary_tracks,
         lst_summary_hratios,
         lst_gene_tracks,
         lst_gene_hratios,
         summary_fsize_x,
         summary_fsize_y,
         gene_fsize_x,
         gene_fsize_y,
         summary_alpha,
         dist_thr, 
         genes,
         n_genes,
         volcano_top_n,
         volcano_fsize_x,
         volcano_fsize_y,
         volcano_subplots_fsize_x,
         volcano_subplots_fsize_y,
         log_odds_fsize_x,
         log_odds_fsize_y,
         output_csv,
         output_all_pos,
         verbose):
    """"Generate summary and individual gene plots for a quick interpretation of the 3D-clustering analysis."""

    startup_message(__version__, "Starting plot generation..")
    logger.info(f"O3D genes-result: {gene_result_path}")
    logger.info(f"O3D positions-result: {pos_result_path}")
    logger.info(f"O3D input mutations: {maf_path}")
    logger.info(f"O3D missense mut prob: {miss_prob_path}")
    logger.info(f"O3D missense seq df: {seq_df_path}")
    logger.info(f"O3D datasets: {datasets_dir}")
    logger.info(f"O3D annotations: {annotations_dir}")
    logger.info(f"Output: {output_dir}")
    logger.info(f"Cohort: {cohort}")
    logger.info(f"Title: {title}")
    logger.info(f"Input mutations including non-missense: {maf_for_nonmiss_path}")
    logger.info(f"Custom summary plot tracks: {lst_summary_tracks}")
    logger.info(f"Custom summary plot h-ratios: {lst_summary_hratios}")
    logger.info(f"Custom gene plots tracks: {lst_gene_tracks}")
    logger.info(f"Custom gene plots h-ratios: {lst_gene_hratios}")
    logger.info(f"Summary plot fsize_x (to be dynamically increased): {summary_fsize_x}")
    logger.info(f"Summary plot fsize_y: {summary_fsize_y}")
    logger.info(f"Gene plots fsize_x: {gene_fsize_x}")
    logger.info(f"Gene plots fsize_y: {gene_fsize_y}")
    logger.info(f"Volcano plot fsize_x: {volcano_fsize_x}")
    logger.info(f"Volcano plot fsize_y: {volcano_fsize_y}")
    logger.info(f"Volcano subplot fsize_x: {volcano_subplots_fsize_x}")
    logger.info(f"Volcano subplot fsize_y: {volcano_subplots_fsize_y}")
    logger.info(f"Log odds plot fsize_x: {log_odds_fsize_x}")
    logger.info(f"Log odds plot fsize_y: {log_odds_fsize_y}")
    logger.info(f"Summary plot score alpha: {summary_alpha}")
    logger.info(f"Threshold for clashing feat: {dist_thr}")
    logger.info(f"Subset of genes: {genes}")
    logger.info(f"Max number of genes to plot: {n_genes}")
    logger.info(f"Volcano plot top associations to annotate: {volcano_top_n}")
    logger.info(f"Output csv file: {bool(output_csv)}")
    logger.info(f"Include non-mutated positions in csv file: {bool(output_all_pos)}")
    logger.info(f"Verbose: {bool(verbose)}")
    logger.info(f'Log path: {os.path.join(output_dir, "log")}')
    logger.info("")

    lst_summary_tracks = parse_lst_tracks(lst_summary_tracks, plot_type="summary")
    lst_gene_tracks = parse_lst_tracks(lst_gene_tracks, plot_type="gene")
    plot_pars = init_plot_pars(summary_fsize_x=summary_fsize_x,
                               summary_fsize_y=summary_fsize_y,
                               gene_fsize_x=gene_fsize_x, 
                               gene_fsize_y=gene_fsize_y, 
                               volcano_fsize_x=volcano_fsize_x,
                               volcano_fsize_y=volcano_fsize_y,
                               volcano_subplots_fsize_x=volcano_subplots_fsize_x,
                               volcano_subplots_fsize_y=volcano_subplots_fsize_y,
                               log_odds_fsize_x=log_odds_fsize_x,
                               log_odds_fsize_y=log_odds_fsize_y,
                               dist_thr=dist_thr,
                               summary_alpha=summary_alpha,
                               lst_summary_tracks=lst_summary_tracks,
                               lst_summary_hratios=lst_summary_hratios,
                               lst_gene_annot=lst_gene_tracks, 
                               lst_gene_hratios=lst_gene_hratios,
                               volcano_top_n=volcano_top_n)

    generate_plots(gene_result_path=gene_result_path,
                  pos_result_path=pos_result_path,
                  maf_path=maf_path,
                  miss_prob_path=miss_prob_path,
                  seq_df_path=seq_df_path,
                  cohort=cohort,
                  datasets_dir=datasets_dir,
                  annotations_dir=annotations_dir,
                  output_dir=output_dir,
                  plot_pars=plot_pars,
                  maf_path_for_nonmiss=maf_for_nonmiss_path,
                  n_genes=n_genes,
                  lst_genes=genes,
                  save_plot=True,
                  show_plot=False,
                  save_csv=output_csv,
                  include_all_pos=output_all_pos,
                  title=title)



# =============================================================================
#                              COMPARATIVE PLOTS
# =============================================================================

# Example:
# oncodrive3D comparative-plot -i /workspace/projects/clustering_3d/o3d_analysys/datasets/output/normal/o3d_output/2024/ -c all_samples_raw -I /workspace/projects/clustering_3d/o3d_analysys/datasets/output/cancer_202404/o3d_output/human_raw/run_2024-05-07_17-46-44 -C TCGA_WXS_BRCA -d /workspace/nobackup/scratch/oncodrive3d/datasets_240506/ -a /workspace/nobackup/scratch/oncodrive3d/annotations_240506/ -o /workspace/projects/clustering_3d/dev_testing/result/plots -v

# TO DO:
# - Test plots with non-missense mutations
# - Enable not morrirs for non-missensem mutations

@oncodrive3D.command(context_settings=dict(help_option_names=['-h', '--help']),
               help="Generate plots to compare two runs of 3D-clustering analysis.")
@click.option("-i", "--o3d_result_dir_1", type=click.Path(exists=True), required=True,
              help="Path to result A directory (including gene- and pos-level result for run A)")
@click.option("-I", "--o3d_result_dir_2", type=click.Path(exists=True), required=True,
              help="Path to result B directory (including gene- and pos-level result for run B)")
@click.option("-c", "--cohort_1", 
              help="Cohort A name (it should match the basename of the files in o3d_result_dir)", type=str, required=True)
@click.option("-C", "--cohort_2", 
              help="Cohort B name (it should match the basename of the files in o3d_result_dir)", type=str, required=True)
@click.option("-d", "--datasets_dir", type=click.Path(exists=True), required=True,
              help="Path to datasets directory")
@click.option("-a", "--annotations_dir", type=click.Path(exists=True), required=True, 
              help="Path to annotations directory")
@click.option("-o", "--output_dir", required=True,
              help="Path to output directory where to save plots")
@click.option("--maf_path_nonmiss_1", type=click.Path(exists=True), 
              help="Path to input mutations file A including non-missense mutations")
@click.option("--maf_path_nonmiss_2", type=click.Path(exists=True), 
              help="Path to input mutations file B including non-missense mutations (e.g., miss_count,miss_prob,score)")
@click.option("--lst_tracks", type=str,
              help="List of tracks to be included in the plots", 
              default="miss_count,miss_prob,score,clusters,ddg,disorder,pacc,ptm,site,sse,pfam,prosite,membrane,motif")
@click.option("--lst_hratios", type=str,
              help="List of float to define horizontal ratio of each track of the plot") 
@click.option("--fsize_x", help="Figure size x-axis", type=float, default=24)
@click.option("--fsize_y", help="Figure size y-axis", type=float, default=12)
@click.option("--dist_thr", help="Threshold of ratios to avoid clashing feature names (e.g., domains and motifs)", type=float, default=0.1)
@click.option("--genes", help="List of genes to be analysed in the report (e.g., --genes TP53,KRAS,PIK3CA)", type=str)
@click.option("--n_genes", help="Top number of genes to be included in the plots", type=int, default=30)
@click.option("--count_mirror", help="Missense mutation count track as mirror image", is_flag=True)
@click.option("--prob_mirror", help="Missense mutation prob track as mirror image", is_flag=True)
@click.option("--score_mirror", help="Clustering score track as mirror image", is_flag=True)
@click.option("-v", "--verbose", help="Verbose", is_flag=True)
@setup_logging_decorator
def comparative_plot(o3d_result_dir_1,
                     cohort_1,
                     o3d_result_dir_2,
                     cohort_2,
                     datasets_dir, 
                     annotations_dir,
                     output_dir,
                     maf_path_nonmiss_1,
                     maf_path_nonmiss_2,
                     lst_tracks,
                     lst_hratios,
                     fsize_x,
                     fsize_y,
                     dist_thr,
                     count_mirror,
                     prob_mirror,
                     score_mirror,
                     genes,
                     n_genes,
                     verbose):
    """"Generate genes comparative plots to comprare two runs of 3D-clustering analysis."""

    startup_message(__version__, "Starting plot generation..")
    logger.info(f"Oncodrive3D result A: {o3d_result_dir_1}")
    logger.info(f"Cohort B: {cohort_1}")
    logger.info(f"Oncodrive3D result B: {o3d_result_dir_2}")
    logger.info(f"Cohort B: {cohort_2}")
    logger.info(f"Oncodrive3D datasets: {datasets_dir}")
    logger.info(f"Oncodrive3D annotations: {annotations_dir}")
    logger.info(f"Output: {output_dir}")
    logger.info(f"Input mutations including non-missense A: {maf_path_nonmiss_1}")
    logger.info(f"Input mutations including non-missense B: {maf_path_nonmiss_2}")
    logger.info(f"Custom gene plots tracks: {lst_tracks}")
    logger.info(f"Custom gene plots h-ratios: {lst_hratios}")
    logger.info(f"Summary plot fsize_x: {fsize_x}")
    logger.info(f"Summary plot fsize_y: {fsize_y}")
    logger.info(f"Threshold for clashing feat: {dist_thr}")
    logger.info(f"Missense mut as mirror image: {bool(count_mirror)}")
    logger.info(f"Missense mut prob as mirror image: {bool(prob_mirror)}")
    logger.info(f"Score as mirror image: {bool(score_mirror)}")
    logger.info(f"Verbose: {bool(verbose)}")
    logger.info(f"Subset of genes: {genes}")
    logger.info(f"Max number of genes to plot: {n_genes}")
    logger.info(f"Verbose: {bool(verbose)}")
    logger.info(f'Log path: {os.path.join(output_dir, "log")}')
    logger.info("")
    
    lst_tracks = parse_lst_tracks(lst_tracks, plot_type="gene")
    plot_pars = init_comp_plot_pars(fsize_x=fsize_x,
                                    fsize_y=fsize_y,
                                    dist_thr=dist_thr,
                                    lst_tracks=lst_tracks,
                                    lst_hratios=lst_hratios,
                                    count_mirror=count_mirror, 
                                    score_mirror=score_mirror,
                                    prob_mirror=prob_mirror)

    generate_comparative_plots(o3d_result_dir_1=o3d_result_dir_1,
                               cohort_1=cohort_1,
                               o3d_result_dir_2=o3d_result_dir_2,
                               cohort_2=cohort_2,
                               datasets_dir=datasets_dir, 
                               annotations_dir=annotations_dir,
                               output_dir=output_dir,
                               plot_pars=plot_pars,
                               maf_path_nonmiss_1=maf_path_nonmiss_1,
                               maf_path_nonmiss_2=maf_path_nonmiss_2,
                               n_genes=n_genes, 
                               lst_genes=genes)

if __name__ == "__main__":
    oncodrive3D()