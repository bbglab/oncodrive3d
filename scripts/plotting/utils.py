import logging
import os

import daiquiri
import subprocess
import click
import sys
import numpy as np
import pandas as pd
import json

from scripts import __logger_name__

logger = daiquiri.getLogger(__logger_name__ + ".annotations.utils")

logging.getLogger('urllib3.connectionpool').setLevel(logging.WARNING)



def get_species(species):
    """
    Simply change species name to accepted format.
    """
    
    if species.capitalize() == "Human" or species.capitalize() == "Homo sapiens":
        species = "Homo sapiens"
    elif species.capitalize() == "Mouse" or species.capitalize() == "Mus musculus": 
        species = "Mus musculus"
    else:
        raise RuntimeError(f"Failed to recognize '{species}' as species. Currently accepted ones are 'Homo sapiens' and 'Mus musculus'. Exiting...")

    return species


def clean_annotations_dir(path: str, loc: str) -> None:
    """
    Clean the annotations directory by removing specific files 
    and subdirectories.

    Args:
        path (str): Path to the directory to be cleaned.
    """

    if loc == "d":

        clean_files = f"rm -rf {os.path.join(path, '*.csv')} {os.path.join(path, '*.tsv')} {os.path.join(path, '*.json')} {os.path.join(path, '.*.txt')}"
        clean_ddg = ["rm", "-rf", os.path.join(path, "stability_change")]
        clean_pdbtool = ["rm", "-rf", os.path.join(path, "pdb_tool")]
        #clean_log = ["rm", "-rf", os.path.join(path, "log")]

        logger.debug(clean_files)
        subprocess.run(clean_files, shell=True)

        logger.debug(' '.join(clean_ddg))
        subprocess.run(clean_ddg)
        
        logger.debug(' '.join(clean_pdbtool))
        subprocess.run(clean_pdbtool)

        # logger.debug(' '.join(clean_log))
        # subprocess.run(clean_log)

    elif loc == "r":
        # TODO: implement cleaning function for output
        pass


def clean_annot_dir(path: str, loc: str = 'd') -> None:
    """
    Clean it upon request if it already exists.

    Args:
        path (str): Path to the directory to be created or cleaned.
    """

    if os.listdir(path) != ['log']:
        logger.warning(f"Directory {path} already exists and is not empty.")

        overwrite = "y" if click.get_current_context().params['yes'] else input("Clean existing directory? (y/n): ")
        while overwrite.lower() not in ["y", "yes", "n", "no"]:
            print("Please choose yes or no")
            overwrite = input("Clean existing directory? (y/n): ")

        if overwrite.lower() in ["y", "yes"]:
            clean_annotations_dir(path, loc)
            logger.info(f"Dataset files in {path} have been removed.")
        else:
            logger.warning(f"Dataset files in {path} have not been removed.")
    else:
        pass
    
    
def get_broad_consequence(list_of_annotations):
    """
    Group variants into broader consequence types.
    """
        
    CONSEQUENCES_LIST = [
        'transcript_ablation',
        'splice_acceptor_variant',
        'splice_donor_variant',
        'stop_gained',
        'frameshift_variant',
        'stop_lost',
        'start_lost',
        'transcript_amplification',
        'inframe_insertion',
        'inframe_deletion',
        'missense_variant',
        'protein_altering_variant',
        'splice_region_variant',
        'splice_donor_5th_base_variant',
        'splice_donor_region_variant',
        'splice_polypyrimidine_tract_variant',
        'incomplete_terminal_codon_variant',
        'start_retained_variant',
        'stop_retained_variant',
        'synonymous_variant',
        'coding_sequence_variant',
        'mature_miRNA_variant',
        '5_prime_UTR_variant',
        '3_prime_UTR_variant',
        'non_coding_transcript_exon_variant',
        'intron_variant',
        'NMD_transcript_variant',
        'non_coding_transcript_variant',
        'upstream_gene_variant',
        'downstream_gene_variant',
        'TFBS_ablation',
        'TFBS_amplification',
        'TF_binding_site_variant',
        'regulatory_region_ablation',
        'regulatory_region_amplification',
        'feature_elongation',
        'regulatory_region_variant',
        'feature_truncation',
        'intergenic_variant'
    ]
    
    GROUPING_DICT = {
        'transcript_ablation': 'nonsense',
        'splice_acceptor_variant': 'nonsense',
        'splice_donor_variant': 'nonsense',
        'stop_gained': 'nonsense',
        'frameshift_variant': 'nonsense',
        'stop_lost': 'nonsense',
        'start_lost': 'nonsense',
        'missense_variant': 'missense',
        'inframe_insertion': 'indel',
        'inframe_deletion': 'indel',
        'splice_donor_variant': 'splicing',
        'splice_acceptor_variant': 'splicing',
        'splice_region_variant': 'splicing',
        'splice_donor_5th_base_variant': 'splicing',
        'splice_donor_region_variant': 'splicing',
        'splice_polypyrimidine_tract_variant': 'splicing',
        'synonymous_variant': 'synonymous',
        'incomplete_terminal_codon_variant': 'synonymous',
        'start_retained_variant': 'synonymous',
        'stop_retained_variant': 'synonymous',
        'protein_altering_variant' : 'protein_altering_variant',
        'transcript_amplification' : 'transcript_amplification', 
        'coding_sequence_variant': 'coding_sequence_variant', 
        'mature_miRNA_variant': 'non_coding_exon_region',
        '5_prime_UTR_variant': 'non_coding_exon_region',
        '3_prime_UTR_variant': 'non_coding_exon_region',
        'non_coding_transcript_exon_variant': 'non_coding_exon_region',
        'NMD_transcript_variant': 'non_coding_exon_region',
        'intron_variant': 'intron_variant',
        'non_coding_transcript_variant' : 'non_coding_transcript_variant',
        'upstream_gene_variant': 'non_genic_variant',
        'downstream_gene_variant': 'non_genic_variant',
        'TFBS_ablation': 'non_genic_variant',
        'TFBS_amplification': 'non_genic_variant',
        'TF_binding_site_variant': 'non_genic_variant',
        'regulatory_region_ablation': 'non_genic_variant',
        'regulatory_region_amplification': 'non_genic_variant',
        'feature_elongation': 'non_genic_variant',
        'regulatory_region_variant': 'non_genic_variant',
        'feature_truncation': 'non_genic_variant',
        'intergenic_variant': 'non_genic_variant',
        '-'  : '-'
    }
    
    consequence_rank_dict = { consequence : rank for rank, consequence in enumerate(CONSEQUENCES_LIST) }
    rank_consequence_dict = { rank : consequence for rank, consequence in enumerate(CONSEQUENCES_LIST) }
    
    list_of_single_annotations = []
    list_of_broad_annotations = []
    for x in list_of_annotations:
        all_consequences = x.split(",")
        all_consequences_ranks = map(lambda x: consequence_rank_dict[x], all_consequences)
        single_consequence = rank_consequence_dict[min(all_consequences_ranks)]
        list_of_single_annotations.append(single_consequence)
        if single_consequence in GROUPING_DICT:
            list_of_broad_annotations.append(GROUPING_DICT[single_consequence])
        else:
            list_of_broad_annotations.append(single_consequence)

    return list_of_broad_annotations


def init_plot_pars(summary_fsize_x=0.4,    # It will be moltiplied for the number of genes
                   summary_fsize_y=8,
                   gene_fsize_x=24, 
                   gene_fsize_y=12, 
                   s_lw=0.2, 
                   sse_fill_width=0.43, 
                   dist_thr=0.1, 
                   summary_alpha=0.3,
                   lst_summary_tracks=None,
                   lst_summary_hratios=None,
                   lst_gene_annot=None, 
                   lst_gene_hratios=None):
    """
    Initialize plotting parameters.
    """
    
    plot_pars = {}

    plot_pars["summary_figsize"] = summary_fsize_x, summary_fsize_y
    plot_pars["figsize"] = gene_fsize_x, gene_fsize_y
    plot_pars["s_lw"] = s_lw
    plot_pars["sse_fill_width"] = sse_fill_width
    plot_pars["dist_thr"] = dist_thr
    plot_pars["summary_alpha"] = summary_alpha

    
    # Summary-plot
    #=============
    
    # Default values
    plot_pars["summary_h_ratios"] = {"score"         : 0.3,
                                     "miss_count"    : 0.2,
                                     "res_count"     : 0.2,
                                     "res_ratio"     : 0.2,
                                     "clusters"      : 0.2}
    
    # Custom values
    if not lst_summary_tracks:
        lst_summary_tracks = plot_pars["summary_h_ratios"].keys()
    if lst_summary_hratios:
        plot_pars["summary_h_ratios"] = {lst_summary_tracks[i] : h_ratio for i, h_ratio in enumerate(lst_summary_hratios)}
    else:
        plot_pars["summary_h_ratios"] = {annot : plot_pars["summary_h_ratios"][annot] for annot in lst_summary_tracks}
    plot_pars["summary_h_ratios"] = {k:v/sum(plot_pars["summary_h_ratios"].values()) for k,v in plot_pars["summary_h_ratios"].items()}
        
    
    # Gene-plots
    #===========
    
    # Default values
    plot_pars["h_ratios"] = {"nonmiss_count" : 0.13,
                             "miss_count"    : 0.13,
                             "miss_prob"     : 0.13,
                             "score"         : 0.13,
                             "pae"           : 0.1,
                             "disorder"      : 0.1,
                             "pacc"          : 0.1,
                             "ddg"           : 0.1,
                             "ptm"           : 0.022,
                             "site"          : 0.022,
                             "clusters"      : 0.04,
                             "sse"           : 0.065,
                             "pfam"          : 0.04,
                             "prosite"       : 0.04,
                             "membrane"      : 0.04,
                             "motif"         : 0.04}
    
    plot_pars["color_cnsq"] = {"splicing" : "C2",
                               "missense" : "C5",
                               "synonymous" : "C9",
                               "coding_sequence_variant" : "C1",
                               "nonsense" : "C6",
                               "intron_variant" : "C7",
                               "indel" : "C8",
                               "protein_altering_variant" : "C3"}

    # Custom values
    if not lst_gene_annot:
        lst_gene_annot = plot_pars["h_ratios"].keys()
    if lst_gene_hratios:
        plot_pars["h_ratios"] = {lst_gene_annot[i] : h_ratio for i, h_ratio in enumerate(lst_gene_hratios)}
    else:
        plot_pars["h_ratios"] = {annot : plot_pars["h_ratios"][annot] for annot in lst_gene_annot}
    
    return plot_pars


def init_comp_plot_pars(fsize_x=24, 
                        fsize_y=12, 
                        s_lw=0.2, 
                        sse_fill_width=0.43, 
                        dist_thr=0.1, 
                        lst_tracks=None,
                        lst_hratios=None,
                        count_mirror=False, 
                        score_mirror=False,
                        prob_mirror=False):
    """
    Initialize plotting parameters.
    """
    
    plot_pars = {}
 
    plot_pars["figsize"] = fsize_x, fsize_y
    plot_pars["s_lw"] = s_lw
    plot_pars["sse_fill_width"] = sse_fill_width
    plot_pars["dist_thr"] = dist_thr
    plot_pars["count_mirror"] = count_mirror
    plot_pars["score_mirror"] = score_mirror
    plot_pars["prob_mirror"] = prob_mirror
        
    # Default values
    plot_pars["h_ratios"] = {"nonmiss_count" : 0.13,
                             "miss_count"    : 0.13,
                             "miss_prob"     : 0.13,
                             "score"         : 0.13,
                             "clusters"      : 0.04,
                             "pae"           : 0.1,
                             "disorder"      : 0.1,
                             "pacc"          : 0.1,
                             "ddg"           : 0.1,
                             "ptm"           : 0.022,
                             "site"          : 0.022,
                             "sse"           : 0.065,
                             "pfam"          : 0.04,
                             "prosite"       : 0.04,
                             "membrane"      : 0.04,
                             "motif"         : 0.04}
    
    plot_pars["color_cnsq"] = {"splicing" : "C2",
                               "missense" : "C5",
                               "synonymous" : "C9",
                               "coding_sequence_variant" : "C1",
                               "nonsense" : "C6",
                               "intron_variant" : "C7",
                               "indel" : "C8",
                               "protein_altering_variant" : "C3"}

    # Custom values
    if not lst_tracks:
        lst_tracks = list(plot_pars["h_ratios"].keys())
    
    if not count_mirror and "nonmiss_count" in lst_tracks:    
        ix = lst_tracks.index("nonmiss_count")
        lst_tracks.insert(ix+1, "nonmiss_count_2")
        plot_pars["h_ratios"]["nonmiss_count_2"] = plot_pars["h_ratios"]["nonmiss_count"] 
        
    if not count_mirror and "miss_count" in lst_tracks:    
        ix = lst_tracks.index("miss_count")
        lst_tracks.insert(ix+1, "miss_count_2")
        plot_pars["h_ratios"]["miss_count_2"] = plot_pars["h_ratios"]["miss_count"]
        
    if not score_mirror and "score" in lst_tracks:    
        ix = lst_tracks.index("score")
        lst_tracks.insert(ix+1, "score_2")
        plot_pars["h_ratios"]["score_2"] = plot_pars["h_ratios"]["score"]
        
    if "clusters" in lst_tracks:
        ix = lst_tracks.index("clusters")
        lst_tracks.insert(ix+1, "clusters_2")
        plot_pars["h_ratios"]["clusters_2"] = plot_pars["h_ratios"]["clusters"]
        
    if "ddg" in lst_tracks:
        ix = lst_tracks.index("ddg")
        lst_tracks.insert(ix+1, "ddg_2")
        plot_pars["h_ratios"]["ddg_2"] = plot_pars["h_ratios"]["ddg"]

    if lst_hratios:
        plot_pars["h_ratios"] = {lst_tracks[i] : h_ratio for i, h_ratio in enumerate(lst_hratios)}
    else:
        plot_pars["h_ratios"] = {annot : plot_pars["h_ratios"][annot] for annot in lst_tracks}   
        
    plot_pars["h_ratios"] = {k:v/sum(plot_pars["h_ratios"].values()) for k,v in plot_pars["h_ratios"].items()}
    
    return plot_pars


def load_o3d_result(o3d_result_path, cohort):
    """
    Load all files generated by 3D clustering analysis of Oncodrive3D.
    """
    
    gene_result_path = f"{o3d_result_path}/{cohort}/{cohort}.3d_clustering_genes.csv"
    pos_result_path = f"{o3d_result_path}/{cohort}/{cohort}.3d_clustering_pos.csv"
    maf_path = f"{o3d_result_path}/{cohort}/{cohort}.mutations.processed.tsv"
    miss_prob_dict_path = f"{o3d_result_path}/{cohort}/{cohort}.miss_prob.processed.json"
    gene_result = pd.read_csv(gene_result_path)
    pos_result = pd.read_csv(pos_result_path)
    maf = pd.read_csv(maf_path, sep="\t")
    miss_prob_dict = json.load(open(miss_prob_dict_path))  
    
    return gene_result, pos_result, maf, miss_prob_dict


def subset_genes_and_ids(genes, 
                         uni_ids, 
                         seq_df, 
                         disorder, 
                         pdb_tool, 
                         uniprot_feat):
    """
    Subset each dataframe by keeping only selected genes and proteins IDs.
    """

    seq_df = seq_df.copy()
    disorder = disorder.copy()
    pdb_tool = pdb_tool.copy()
    uniprot_feat = uniprot_feat.copy()
    # Filter genes in the other df
    seq_df = seq_df[seq_df["Gene"].isin(genes)]
    disorder = disorder[disorder["Uniprot_ID"].isin(uni_ids)].reset_index(drop=True) 
    pdb_tool = pdb_tool[pdb_tool["Uniprot_ID"].isin(uni_ids)].reset_index(drop=True)
    uniprot_feat = uniprot_feat[uniprot_feat["Gene"].isin(genes)]
    
    return seq_df, disorder, pdb_tool, uniprot_feat


def filter_o3d_result(gene_result, pos_result, n_genes=None, lst_genes=None):
    """
    Subset gene-level and position-level Oncodrive3D result.
    """

    if isinstance(lst_genes, str):
        lst_genes = lst_genes.replace(" ", "")
        lst_genes = lst_genes.split(",")
        gene_result = gene_result[gene_result["Gene"].isin(lst_genes)]    
    gene_result = gene_result[gene_result["Status"] == "Processed"]
    if n_genes:
        gene_result = gene_result[:n_genes]                                     
    uni_ids = gene_result.Uniprot_ID.values
    genes = gene_result.Gene.values   
    pos_result = pos_result[pos_result["Gene"].isin(genes)]   
    
    return gene_result, pos_result, genes, uni_ids


def get_enriched_result(pos_result_gene, 
                        disorder_gene, 
                        pdb_tool_gene, 
                        seq_df):
    """
    Add annotations to Oncodrive3D result to return an annotated tsv. 
    """

    pos_result_gene = pos_result_gene.copy()
    # DDG
    pos_result_gene.loc[pos_result_gene["Mut_in_res"] == 0, "DDG"] = np.nan
    # Disorder
    pos_result_gene = pos_result_gene.merge(disorder_gene, how="left", on=["Pos"])
    pos_result_gene = pos_result_gene.rename(columns={"Confidence" : "pLDDT_res"})
    # PDB_Tool
    pos_result_gene = pos_result_gene.rename(columns={"AF_F" : "F"}).merge(
            pdb_tool_gene, on=["Res", "Uniprot_ID", "Pos", "F"])
    # Transcript and gene IDs
    pos_result_gene = pos_result_gene.merge(
        seq_df[["Gene", "Uniprot_ID", "Ens_Gene_ID", "Ens_Transcr_ID"]], 
        how="left", on=["Uniprot_ID"])

    return pos_result_gene
    
    
def save_annotated_pos_result(pos_result, 
                               annot_pos_result, 
                               pfam_processed, 
                               output_dir, 
                               run_name, 
                               output_all_pos=False):
    """
    Save the annotated pos-level result.
    """
    
    # Do not include non-mutated positions (default)
    if output_all_pos == False:
        annot_pos_result = annot_pos_result[annot_pos_result["Mut_in_res"] != 0].reset_index(drop=True)
        
    # Merge with 'original' one to retrieve dropped cols
    output_pos_result = os.path.join(output_dir, f"{run_name}.3d_clustering_pos.annotated.tsv")
    output_pfam = os.path.join(output_dir, f"{run_name}.pfam.tsv")
    annot_pos_result = pos_result.drop(columns=["F", "pLDDT_res"]).merge(
        annot_pos_result[["Gene", 
                          "Uniprot_ID",
                          "F", 
                          "Ens_Gene_ID", 
                          "Ens_Transcr_ID", 
                          "Pos",
                          "Res", 
                          "pLDDT_res", 
                          "SSE", 
                          "pACC", 
                          "DDG"]],
        how="right", on=["Gene", "Uniprot_ID", "Pos"])
    annot_pos_result = annot_pos_result.sort_values(["Gene", "Pos"])
    
    # Fill the NA of the non-mutated positions in features
    if annot_pos_result["Cancer"].isnull().all():
        annot_pos_result["Cancer"] = np.nan
    else:
        annot_pos_result["Cancer"] = annot_pos_result["Cancer"].dropna().unique()[0]
    if annot_pos_result["Cohort"].isnull().all():
        annot_pos_result["Cohort"] = np.nan
    else:
        annot_pos_result["Cohort"] = annot_pos_result["Cohort"].dropna().unique()[0]
    annot_pos_result["Mut_in_res"] = annot_pos_result["Mut_in_res"].fillna(0)
    for gene in annot_pos_result.Gene.unique():
        mut_in_gene = annot_pos_result.loc[annot_pos_result["Gene"] == gene, "Mut_in_gene"].dropna().unique()[0]
        tot_samples = annot_pos_result.loc[annot_pos_result["Gene"] == gene, "Tot_samples"].dropna().unique()
        if tot_samples: 
            tot_samples = tot_samples[0]
        else:
            tot_samples = np.nan
        annot_pos_result.loc[annot_pos_result["Gene"] == gene, "Mut_in_gene"] = mut_in_gene
        annot_pos_result.loc[annot_pos_result["Gene"] == gene, "Tot_samples"] = tot_samples

    # Save
    annot_pos_result.to_csv(output_pos_result, sep="\t", index=False)
    pfam_processed.to_csv(output_pfam, sep="\t", index=False)
    logger.info(f"Saved annotated position-level result to {output_pos_result}")
    logger.info(f"Saved Pfam annotations to {output_pfam}")
    
    
def parse_lst_tracks(lst, plot_type):
    """
    Parse the list of tracks from click arg.
    """
    
    summary_tracks = ["score",
                      "miss_count",
                      "res_count",
                      "res_ratio",
                      "clusters"]
    
    gene_tracks = ["nonmiss_count",
                  "miss_count",
                  "miss_prob",
                  "score",
                  "pae",
                  "disorder",
                  "pacc",
                  "ddg",
                  "ptm",
                  "site",
                  "clusters",
                  "sse",
                  "pfam",
                  "prosite",
                  "membrane",
                  "motif"]

    if plot_type == "summary":
        available_tracks = summary_tracks
    elif plot_type == "gene":
        available_tracks = gene_tracks
    lst = lst.split(",")
    
    is_valid = np.array([track not in available_tracks for track in lst])
    if is_valid.any():
        invalid_tracks = list(np.array(lst)[np.where(is_valid)])
        logger.error(f"One or more track names for {plot_type} plot are not accepted: {invalid_tracks}")
        logger.error(f"Available track names are: {available_tracks}")
        logger.error(f"Exiting..")
        sys.exit(1)
        
    return lst