import logging
import os

import daiquiri
import subprocess
import click
from scipy import interpolate
import numpy as np

from scripts import __logger_name__

logger = daiquiri.getLogger(__logger_name__ + ".annotations.utils")

logging.getLogger('urllib3.connectionpool').setLevel(logging.WARNING)



# Utils
# =====


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


def init_annotations(annotations):
    """
    Save string of annotations provided as argument as dictionary of bolean.
    """
    
    plot_annot = {}
    plot_annot["nonmiss_count"] = False
    plot_annot["pae"] = False
    plot_annot["disorder"] = False
    plot_annot["pacc"] = False
    plot_annot["ddg"] = False
    plot_annot["ptm"] = False
    plot_annot["site"] = False
    plot_annot["clusters"] = False
    plot_annot["sse"] = False
    plot_annot["pfam"] = False   
    plot_annot["prosite"] = False  
    plot_annot["membrane"] = False
    plot_annot["motif"] = False
    
    if isinstance(annotations, str):
        lst_annotations = [annot.lower() for annot in annotations.replace(" ", "").split(",")]
        if "all" in lst_annotations:
            for annot in plot_annot.keys():
                plot_annot[annot] = True
        elif "none" in lst_annotations:
            for annot in plot_annot.keys():
                plot_annot[annot] = False
        else:
            for annot in lst_annotations:
                if annot in plot_annot:
                    plot_annot[annot] = True
                else:
                    print(f"{annot} is not recognized as valid annotation.")  
        lst_annotations = [k for k,v in plot_annot.items() if v == True]
        if len(lst_annotations) > 0:
            logger.info(f"The following annotations will be included: {lst_annotations}")
        else:
            logger.info("No annotations will be included.")
    else:
        logger.info("No extra annotations will be included in the plots.")
    
    return plot_annot


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


def filter_o3d_result(gene_result, pos_result, n_genes, lst_genes, non_significant):
    """
    Subset gene-level and position-level Oncodrive3D result.
    """

    if isinstance(lst_genes, str):
        lst_genes = lst_genes.replace(" ", "")
        lst_genes = lst_genes.split(",")
        gene_result = gene_result[[gene in lst_genes for gene in gene_result["Gene"].values]]    
    if non_significant == False:
        gene_result = gene_result[gene_result["C_gene"] == 1]
    gene_result = gene_result[gene_result["Status"] == "Processed"]
    gene_result = gene_result[:n_genes]                                     
    uni_ids = gene_result.Uniprot_ID.values
    genes = gene_result.Gene.values   
    pos_result = pos_result[[gene in genes for gene in pos_result["Gene"].values]]   
    
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


def interpolate_x_y(x, y):
    """
    Interpolate x, y to get a completely filled plot without gaps.
    """

    step = (0.01 * len(x)) / 393    # Step coefficient tuned on TP53

    # Do nothing if the step is larger than the original ones (positions)
    if step >= 1:
        
        return x, y

    # Interpolate between positions and values
    else:
        x = x.copy()
        y = y.copy()
        x = np.array(x)
        y = np.array(y)
        
        # Densify data for filled color plot
        f = interpolate.interp1d(x, y)
        x = np.arange(x[0], x[-1], 0.01)
        y = f(x)
    
        return x, y 
    
    
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
        tot_samples = annot_pos_result.loc[annot_pos_result["Gene"] == gene, "Tot_samples"].dropna().unique()[0]
        annot_pos_result.loc[annot_pos_result["Gene"] == gene, "Mut_in_gene"] = mut_in_gene
        annot_pos_result.loc[annot_pos_result["Gene"] == gene, "Tot_samples"] = tot_samples

    # Save
    annot_pos_result.to_csv(output_pos_result, sep="\t", index=False)
    pfam_processed.to_csv(output_pfam, sep="\t", index=False)
    logger.info(f"Saved annotated position-level result to {output_pos_result}")
    logger.info(f"Saved Pfam annotations to {output_pfam}")