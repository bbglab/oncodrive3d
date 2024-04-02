import pandas as pd

import daiquiri
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import os
import json
import colorcet as cc
from scripts.run.utils import parse_maf_input
from scripts.plotting.utils import get_broad_consequence, save_annotated_pos_result, interpolate_x_y
from scripts.plotting.utils import get_enriched_result, filter_o3d_result, subset_genes_and_ids
from scripts.run.mutability import init_mutabilities_module
from scripts.run.miss_mut_prob import get_miss_mut_prob_dict, mut_rate_vec_to_dict, get_unif_gene_miss_prob

from scripts import __logger_name__

logger = daiquiri.getLogger(__logger_name__ + ".annotations.plot")


# Missense mutations probability
# ==============================

def get_miss_mut_prob_for_plot(mut_profile_path, mutability_config_path, seq_df):

    if mutability_config_path is not None:
        # Compute dict from mutability
        logger.debug(f"Computing missense mut probabilities using mutabilities...")
        mutab_config = json.load(open(mutability_config_path))
        init_mutabilities_module(mutab_config)
        seq_df = seq_df[seq_df["Reference_info"] == 1]   
        seq_df['Exons_coord'] = seq_df['Exons_coord'].apply(eval)  
        miss_prob_dict = get_miss_mut_prob_dict(mut_rate_dict=None, 
                                                seq_df=seq_df,
                                                mutability=True, 
                                                mutability_config=mutab_config)
    elif mut_profile_path is not None:
        # Compute dict from mut profile of the cohort and dna sequences
        mut_profile = json.load(open(mut_profile_path))
        logger.debug(f"Computing missense mut probabilities...")
        if not isinstance(mut_profile, dict):
            mut_profile = mut_rate_vec_to_dict(mut_profile)
        miss_prob_dict = get_miss_mut_prob_dict(mut_rate_dict=mut_profile, 
                                                seq_df=seq_df)
    else:
        logger.debug(f"Mutation profile not provided: Uniform distribution will be used for scoring and simulations.")
        miss_prob_dict = None

    return miss_prob_dict


# Summary plots
# =============

def create_plot_dir(directory_path):
    
    # Check if the directory already exists
    if not os.path.exists(directory_path):
        # If not, create the directory
        os.makedirs(directory_path)
        logger.debug(f"Directory '{directory_path}' created successfully.")
    else:
        logger.debug(f"Directory '{directory_path}' already exists.")
        
        
def get_summary_counts(gene_result, pos_result, seq_df):
    """
    Get dataframes including the counts required to generate the summary plots.
    """

    # Df with mut count
    count_mut_gene_hit = gene_result[["Gene", "Clust_mut"]].rename(columns={"Clust_mut" : "Mut_in_gene"})
    count_mut_gene_hit["C"] = "In cluster"
    count_mut_gene_not = pd.DataFrame({"Gene" : count_mut_gene_hit["Gene"].values,
                                    "Mut_in_gene" : gene_result.apply(lambda x: x["Mut_in_gene"] - x["Clust_mut"], axis=1)})
    count_mut_gene_not["C"] = "Not in cluster"
    count_mut_gene_df = pd.concat((count_mut_gene_hit, count_mut_gene_not)).sort_values("Gene").rename(columns={"Mut_in_gene" : "Count"})
    count_mut_gene_df = count_mut_gene_df.sort_values(["C", "Count"], ascending=False).reset_index(drop=True)

    # Df with pos count
    pos_result_not = pos_result[pos_result["C"] == 0]
    if len(pos_result_not) > 0: 
        pos_result_not = pos_result_not.groupby("Gene").apply(len)
        pos_result_not = pos_result_not.reset_index().rename(columns={0 : "Count"})
        pos_result_not["C"] = "Not significant"
    else:
        pos_result_not = pd.DataFrame(columns=["Gene", "Count", "C"])
    pos_result_hit = pos_result[pos_result["C"] == 1]
    if len(pos_result_hit) > 0:   
        pos_result_hit = pos_result_hit.groupby("Gene").apply(len)
        pos_result_hit = pos_result_hit.reset_index().rename(columns={0 : "Count"})
        pos_result_hit["C"] = "Significant"
    else:
        pos_result_hit = pd.DataFrame(columns=["Gene", "Count", "C"])

    count_pos_df = pd.concat((pos_result_hit, pos_result_not)).sort_values("Gene")
    count_pos_df = count_pos_df.sort_values("C", ascending=False).reset_index(drop=True)

    # Df with cluster count
    cluster_df = pos_result.groupby("Gene").max("Cluster").Cluster.reset_index()
    cluster_df["Cluster"] = cluster_df["Cluster"] + 1

    # Df with size
    size_df = pd.DataFrame(seq_df.apply(lambda x: (x.Gene, len(x.Seq)), axis=1).to_list())
    size_df.columns = "Gene", "Size"
    
    return count_mut_gene_df, count_pos_df, cluster_df, size_df
        

def summary_plot(gene_result, 
                 pos_result, 
                 count_mut_gene_df, 
                 count_pos_df, 
                 cluster_df,
                 size_df,
                 output_dir,
                 cohort):        

    # Plot
    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, 
                                             figsize=(10, 12), 
                                             sharex=True, 
                                             gridspec_kw={'hspace': 0.05, 
                                                          'height_ratios': [0.2, 0.2, 0.2, 0.2, 0.3]})

    sns.barplot(x='Gene', y='Size', data=size_df, order=gene_result.Gene, ax=ax1, color="lightgray", ec="black", lw=0.5)
    sns.barplot(x='Gene', y='Count', data=count_mut_gene_df, order=gene_result.Gene, hue="C", ax=ax2, ec="black", lw=0.5)
    sns.barplot(x='Gene', y='Count', data=count_pos_df, order=gene_result.Gene, hue="C", ax=ax3,
                palette=sns.color_palette("tab10", n_colors=2), ec="black", lw=0.5)
    sns.barplot(x='Gene', y='Cluster', data=cluster_df, order=gene_result.Gene, ax=ax4, color="lightgray", ec="black", lw=0.5)

    pos_result = pos_result.copy()
    pos_result["C"] = pos_result.C.map({1 : "Significant", 0 : "Not significant"})
    sns.boxplot(x='Gene', y='Ratio_obs_sim', data=pos_result, order=gene_result.Gene, color="lightgray", showfliers=False, ax=ax5)
    sns.stripplot(x='Gene', y='Ratio_obs_sim', data=pos_result, hue="C" ,jitter=True, size=6, alpha=0.5, order=gene_result.Gene.values, 
                palette=sns.color_palette("tab10", n_colors=2), ax=ax5)

    # Add "X" or "O" for significant and not
    for i, gene in enumerate(gene_result['Gene']):
        median_value = pos_result[pos_result['Gene'] == gene]['Ratio_obs_sim'].median()
        is_significant = gene_result[gene_result["Gene"] == gene].C_gene.values[0]
        text_position = median_value - 0.05
        text_symbol = "*" if is_significant else ""
        ax5.text(i, text_position, text_symbol, ha='center', va='center', fontsize=22, fontweight='bold', color="Black")

    # Details
    fig.suptitle(f'Oncodrive3D output', fontsize=16)
    ax1.set_xlabel(None)
    ax2.set_xlabel(None)
    ax3.set_xlabel(None)
    ax4.set_xlabel(None)
    ax1.set_ylabel('Protein length', fontsize=14)
    ax2.set_ylabel('Missense\nmutations #', fontsize=14)
    ax3.set_ylabel('Mutated\nresidues #', fontsize=14)
    ax4.set_ylabel('Clusters #', fontsize=14)
    ax5.set_ylabel('O3D score\n(obs / μ simulated)', fontsize=14)
    ax2.legend(fontsize=12, loc="upper right")
    ax3.legend(fontsize=12, loc="upper right")
    ax5.legend(fontsize=12)
    plt.xticks(rotation=45, rotation_mode="anchor", ha='right', fontsize=12)
    plt.subplots_adjust(top=0.95) 
    
    # Save
    filename = f"{cohort}.summary_plot.png"
    output_path = os.path.join(output_dir, filename)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    logger.debug(f"Saved {output_path}")
    plt.close()

    
# Gene plots
# ==========

def capitalize(string):
    words = string.split("_")
    words[0] = words[0].capitalize()

    return ' '.join(words)

                           
def get_nonmiss_mut(path_to_maf):
    """
    Get non missense mutations from MAF file.
    """
    try:
        maf_nonmiss = pd.read_csv(path_to_maf, sep="\t", dtype={'Chromosome': str})
        maf_nonmiss = maf_nonmiss[maf_nonmiss["Protein_position"] != "-"]                                            ## TODO: Fix it for alternative MAF (see cancer) 
        maf_nonmiss = maf_nonmiss[~(maf_nonmiss['Consequence'].str.contains('Missense_Mutation')
                                    | maf_nonmiss['Consequence'].str.contains('missense_variant'))]
        maf_nonmiss = maf_nonmiss[["SYMBOL", 
                                "Consequence", 
                                "Protein_position"]].rename(
                                    columns={"SYMBOL" : "Gene", 
                                                "Protein_position" : "Pos"}).reset_index(drop=True)
                                
        # Parse the consequence with multiple elements and get broader categories
        maf_nonmiss["Consequence"] = get_broad_consequence(maf_nonmiss["Consequence"])
                
        return maf_nonmiss
    
    except Exception as e:
        logger.warning("Can't parse non-missense mutation from MAF file: The track will not be included...")           ## TODO: Capture the error better


def avg_per_pos_ddg(pos_result_gene, ddg_prot, maf_gene):
    """
    Compute per-position average stability change upon mutations (DDG).
    """
    
    ddg_vec = np.repeat(0., len(pos_result_gene))
    for pos, group in maf_gene.groupby('Pos'):
        pos = str(pos)
        obs_mut = group.Mut
        if pos in ddg_prot:
            ddg_pos = ddg_prot[pos]
            ddg_pos = np.mean([ddg_pos[mut] for mut in obs_mut])
            ddg_vec[int(pos)-1] = ddg_pos

    return ddg_vec


def check_near_domains(uni_feat_gene, pfam=True, dist_thr=0.05):
    """
    Check if two domains could be closer to each other 
    than allowed threshold (ratio of protein size).
    """
        
    near_domains = False
    uni_feat_gene = uni_feat_gene.copy()
    uni_feat_gene = uni_feat_gene[uni_feat_gene["Type"] == "DOMAIN"]
    if pfam:
        uni_feat_gene = uni_feat_gene[uni_feat_gene["Evidence"] == "Pfam"]
    else:
        uni_feat_gene = uni_feat_gene[uni_feat_gene["Evidence"] != "Pfam"]
    uni_feat_gene = uni_feat_gene.drop_duplicates(subset='Description', keep='first')   ## Does it need it?
    mid_pos = (uni_feat_gene.Begin + uni_feat_gene.End) / 2
    mid_pos_norm = (mid_pos / mid_pos.max()).values
    
    for i in range(len(mid_pos_norm)):
        for j in range(i + 1, len(mid_pos_norm)):
            diff = abs(mid_pos_norm[i] - mid_pos_norm[j])
            if diff < dist_thr:
                near_domains = True

    return near_domains


def check_near_motifs(uni_feat_gene, dist_thr=0.05):
    """
    Check if two motifs could be closer to each other 
    than allowed threshold (ratio of protein size).
    """
        
    near_motif = False
    uni_feat_gene = uni_feat_gene.copy()
    uni_feat_gene = uni_feat_gene[uni_feat_gene["Type"] == "MOTIF"]
    uni_feat_gene = uni_feat_gene.drop_duplicates(subset='Full_description', keep='first')   ## Does it need it?
    mid_pos = (uni_feat_gene.Begin + uni_feat_gene.End) / 2
    mid_pos_norm = (mid_pos / mid_pos.max()).values
    
    for i in range(len(mid_pos_norm)):
        for j in range(i + 1, len(mid_pos_norm)):
            diff = abs(mid_pos_norm[i] - mid_pos_norm[j])
            if diff < dist_thr:
                near_motif = True

    return near_motif


def parse_pos_result_for_genes_plot(pos_result_gene):
    """
    Get mut count and score divided by Oncodriv3D 
    result: significant, not significant, significant extended
    (mutation in a non-significant residue that contribute to 
    mutations in the volume of a significant one/s).
    """
    
    pos_result_gene = pos_result_gene.copy()
    pos_result_gene = pos_result_gene[["Pos", "Mut_in_res", "Mut_in_vol", "Ratio_obs_sim", "C", "C_ext", "pval", "Cluster", "PAE_vol"]]
    pos_result_gene["C"] = pos_result_gene.apply(
        lambda x: 1 if (x["C"] == 1) & (x["C_ext"] == 0) else 2 if (x["C"] == 1) & (x["C_ext"] == 1) else 0, axis=1)

    pos_result_gene_hit = pos_result_gene[pos_result_gene["C"] == 1]
    pos_result_gene_not = pos_result_gene[pos_result_gene["C"] == 0]
    pos_result_gene_ext = pos_result_gene[pos_result_gene["C"] == 2]
    max_mut = np.max(pos_result_gene["Mut_in_res"].values)

    pos_hit = pos_result_gene_hit["Pos"].values
    pos_not = pos_result_gene_not["Pos"].values
    pos_ext = pos_result_gene_ext["Pos"].values
    pos_hit_mut = pos_result_gene_hit["Mut_in_res"].values
    pos_not_mut = pos_result_gene_not["Mut_in_res"].values
    pos_ext_mut = pos_result_gene_ext["Mut_in_res"].values
    pos_hit_score = pos_result_gene_hit["Ratio_obs_sim"].values
    pos_not_score = pos_result_gene_not["Ratio_obs_sim"].values
    pos_ext_score = pos_result_gene_ext["Ratio_obs_sim"].values
    
    return [pos_result_gene, max_mut, 
            pos_hit, pos_not, pos_ext, 
            pos_hit_mut, pos_not_mut, pos_ext_mut, 
            pos_hit_score, pos_not_score, pos_ext_score]
    
    
def get_count_for_genes_plot(maf, maf_nonmiss, gene, non_missense_count=False):
    """
    Get missense and non-missense mutations count.
    """
    
    mut_count = maf.value_counts("Pos").reset_index()
    mut_count = mut_count.rename(columns={0 : "Count"})
    if non_missense_count:
        maf_nonmiss_gene = maf_nonmiss[maf_nonmiss["Gene"] == gene]
        mut_count_nonmiss = maf_nonmiss_gene.groupby("Consequence").value_counts("Pos").reset_index()
        mut_count_nonmiss = mut_count_nonmiss.rename(columns={0 : "Count"})
        # If there is more than one position affected, take the first one
        ix_more_than_one_pos = mut_count_nonmiss.apply(lambda x: len(x["Pos"].split("-")), axis=1) > 1
        mut_count_nonmiss.loc[ix_more_than_one_pos, "Pos"] = mut_count_nonmiss.loc[ix_more_than_one_pos].apply(lambda x: x["Pos"].split("-")[0], axis=1)
        # Filter non-numerical Pos and get count
        mut_count_nonmiss = mut_count_nonmiss[mut_count_nonmiss["Pos"].apply(lambda x: x.isdigit() or x.isnumeric())]
        mut_count_nonmiss["Pos"] = mut_count_nonmiss["Pos"].astype(int)
    else:
        mut_count_nonmiss = None

    return mut_count, mut_count_nonmiss


def get_score_for_genes_plot(pos_result_gene, mut_count, prob_vec):
    """
    Add any non-mutated position to the pos_result df, get 
    per-position score and normalized score.
    """

    pos_result_gene = pos_result_gene.copy()
    score_vec = []
    for pos in range(1, len(prob_vec)+1):

        # Mut count
        if pos in mut_count.Pos.values:
            if pos not in pos_result_gene.Pos.values:
                logger.error("Position in MAF not found in position-level O3D result: Check that MAF and O3D result are matching!")
            score = pos_result_gene.loc[pos_result_gene["Pos"] == pos, "Ratio_obs_sim"].values[0]
        else:
            score = 0
            row_gene = pd.DataFrame({'Pos': [pos], 'Mut_in_res': [0], 'Ratio_obs_sim': [np.nan], 'C': [np.nan]})
            pos_result_gene = pd.concat([pos_result_gene, row_gene])

        score_vec.append(score)

    pos_result_gene = pos_result_gene.sort_values("Pos").reset_index(drop=True)
    
    # Normalize score
    if np.isnan(score_vec).any():
        score_vec = pd.Series(score_vec).fillna(max(score_vec)).values
    score_norm_vec = np.array(score_vec) / sum(score_vec) 
    
    return pos_result_gene, score_vec, score_norm_vec


def get_id_annotations(uni_id, pos_result_gene, maf_gene, annotations_dir, disorder, pdb_tool, uniprot_feat):
    """
    Get the annotations for a specific protein ID.
    """
    
    pos_result_gene = pos_result_gene.copy()
    disorder_gene = disorder[disorder["Uniprot_ID"] == uni_id].reset_index(drop=True)
    pdb_tool_gene = pdb_tool[pdb_tool["Uniprot_ID"] == uni_id].reset_index(drop=True)
    uni_feat_gene = uniprot_feat[uniprot_feat["Uniprot_ID"] == uni_id].reset_index(drop=True)
    ddg_path = os.path.join(annotations_dir, "stability_change", f"{uni_id}_ddg.json")
    if os.path.isfile(ddg_path):
        ddg = json.load(open(ddg_path))
        ddg_vec = avg_per_pos_ddg(pos_result_gene, ddg, maf_gene)
        pos_result_gene["DDG"] = ddg_vec
    else:
        pos_result_gene["DDG"] = np.nan

    # Avoid duplicates (Uniprot IDs mapping to the different gene names)
    uni_feat_gene = uni_feat_gene.drop(columns=["Gene", "Ens_Transcr_ID", "Ens_Gene_ID"]).drop_duplicates()
    
    return pos_result_gene, disorder_gene, pdb_tool_gene, uni_feat_gene


def set_axes_arg(pos_result_gene, plot_pars, plot_annot_gene, uni_feat_gene):
    """
    Adjust the height ratio (and the number) of tracks to include 
    in the plot. 
    """
    
    dx = 0
    h_ratios = plot_pars["h_ratios"].copy()

    if plot_annot_gene["nonmiss_count"] == False:
        del h_ratios[0]
        dx += 1
    if np.isnan(pos_result_gene["PAE_vol"]).all() or plot_annot_gene["pae"] == False:
        del h_ratios[3-dx]
        dx += 1
        plot_annot_gene["pae"] = False
    if plot_annot_gene["disorder"] == False:
        del h_ratios[4-dx]
        dx += 1
    if plot_annot_gene["pacc"] == False:
        del h_ratios[5-dx]
        dx += 1
    if plot_annot_gene["ddg"] == False:
        del h_ratios[6-dx]
        dx += 1

    if plot_annot_gene["ptm"] == False or len(uni_feat_gene[uni_feat_gene["Type"] == "PTM"]) == 0:
        del h_ratios[7-dx]
        dx += 1
        plot_annot_gene["ptm"] = False
    else:
        stracks = len(uni_feat_gene[uni_feat_gene["Type"] == "PTM"].Description.unique())
        h_ratios[7-dx] = h_ratios[7-dx] * stracks
 
    if plot_annot_gene["site"] == False or len(uni_feat_gene[uni_feat_gene["Type"] == "SITE"]) == 0:
        del h_ratios[8-dx]
        dx += 1
        plot_annot_gene["site"] = False
    else:
        stracks = len(uni_feat_gene[uni_feat_gene["Type"] == "SITE"].Description.unique())
        h_ratios[8-dx] = h_ratios[8-dx] * stracks

    
    if plot_annot_gene["clusters"] == False:
        del h_ratios[9-dx]
        dx += 1
    if plot_annot_gene["sse"] == False:
        del h_ratios[10-dx]
        dx += 1

    if plot_annot_gene["pfam"] == False or len(uni_feat_gene[(uni_feat_gene["Type"] == "DOMAIN") & 
                                          (uni_feat_gene["Evidence"] == "Pfam")]) == 0:
        del h_ratios[11-dx]
        dx += 1
        near_pfam = False
        plot_annot_gene["pfam"] = False
    else:
        # Make the track for the Pfam domains larger if domains are too close
        near_pfam = check_near_domains(uni_feat_gene, pfam=True, dist_thr=plot_pars["dist_thr"])
        if near_pfam:
            h_ratios[11-dx] = 0.08
            
    if plot_annot_gene["prosite"] == False or len(uni_feat_gene[(uni_feat_gene["Type"] == "DOMAIN") & 
                                          (uni_feat_gene["Evidence"] != "Pfam")]) == 0:
        del h_ratios[12-dx]
        dx += 1
        near_prosite = False
        plot_annot_gene["prosite"] = False
    else:
        near_prosite = check_near_domains(uni_feat_gene, pfam=False, dist_thr=plot_pars["dist_thr"])
        if near_prosite:
            h_ratios[12-dx] = 0.08
            
    if plot_annot_gene["membrane"] == False or len(uni_feat_gene[uni_feat_gene["Type"] == "MEMBRANE"]) == 0:
        del h_ratios[13-dx]
        dx += 1
        plot_annot_gene["membrane"] = False
    else:
        stracks = len(uni_feat_gene[uni_feat_gene["Type"] == "MEMBRANE"].Description.unique())
        h_ratios[13-dx] = h_ratios[13-dx] * stracks
    
    if plot_annot_gene["motif"] == False or len(uni_feat_gene[uni_feat_gene["Type"] == "MOTIF"]) == 0:
        del h_ratios[14-dx]
        dx += 1
        near_motif = False
        plot_annot_gene["motif"] = False
    else:
        near_motif = check_near_domains(uni_feat_gene, dist_thr=0.1)
        if near_motif:
            h_ratios[14-dx] = 0.08
        
    return h_ratios, plot_annot_gene, near_pfam, near_prosite, near_motif


def genes_plots(gene_result, 
                pos_result, 
                seq_df,
                maf,
                miss_prob_dict,
                output_dir,
                cohort,
                annotations_dir,
                disorder,
                uniprot_feat,
                pdb_tool,
                maf_nonmiss,
                plot_annot,
                plot_pars,
                size_df,
                output_tsv=False):   
    """
    Generate a diagnostic plot for each gene showing Oncodrive3D 
    results and annotated features.
    """
    
    # Parameters
    # ==========
    
    annotated_result_lst = []
    uni_feat_result_lst = []
    for j, gene in enumerate(gene_result["Gene"].values):
      
        # Load and parse
        # ==============
        
        # IDs
        uni_id = seq_df[seq_df["Gene"] == gene].Uniprot_ID.values[0]
        af_f = seq_df[seq_df["Gene"] == gene].F.values[0]
        maf_gene = maf[maf["Gene"] == gene]

        # Parse
        pos_result_gene = pos_result[pos_result["Gene"] == gene].sort_values("Pos").reset_index(drop=True)

        if len(pos_result_gene) > 0:
            pos_result_gene = pos_result_gene[["Pos", "Mut_in_res", "Mut_in_vol", 
                                               "Ratio_obs_sim", "C", "C_ext", 
                                               "pval", "Cluster", "PAE_vol"]]
            [pos_result_gene, max_mut, 
             pos_hit, pos_not, pos_ext, 
             pos_hit_mut, pos_not_mut, pos_ext_mut, 
             pos_hit_score, pos_not_score, pos_ext_score] = parse_pos_result_for_genes_plot(pos_result_gene)
            
            # Counts
            mut_count, mut_count_nonmiss = get_count_for_genes_plot(maf_gene, 
                                                                    maf_nonmiss, 
                                                                    gene, 
                                                                    plot_annot["nonmiss_count"])
            
            # Get prob vec
            prob_vec = miss_prob_dict[f"{uni_id}-F{af_f}"]                         # TODO: If none, use uniform        <-------------------------- TODO
        
            # Get per-pos score and normalize score
            pos_result_gene, score_vec, score_norm_vec = get_score_for_genes_plot(pos_result_gene, 
                                                                                  mut_count, 
                                                                                  prob_vec)

            # Get annotations
            pos_result_gene, disorder_gene, pdb_tool_gene, uni_feat_gene = get_id_annotations(uni_id, 
                                                                                               pos_result_gene, 
                                                                                               maf_gene, 
                                                                                               annotations_dir, 
                                                                                               disorder, 
                                                                                               pdb_tool, 
                                                                                               uniprot_feat)

            # Generate plot
            # ============= 
                
            ax = 0
            gene_len = size_df[size_df["Gene"] == gene].Size.values[0]
            plot_annot_gene = plot_annot.copy()
            if mut_count_nonmiss is None:
                plot_annot_gene["nonmiss_count"] = False
            h_ratios, plot_annot_gene, near_pfam, near_prosite, near_motif = set_axes_arg(pos_result_gene, plot_pars, plot_annot_gene, uni_feat_gene)

            ntracks = len(h_ratios)
            h_ratios = np.array(h_ratios) / sum(h_ratios)

            fig, axes = plt.subplots(ntracks, 1, 
                                    figsize=plot_pars["figsize"], sharex=True, 
                                    gridspec_kw={'hspace': 0.1, 
                                                'height_ratios': h_ratios})
                
            # Plot for Non-missense mut track   
            # -------------------------------
            if plot_annot_gene["nonmiss_count"]:
                try:
                    if len(mut_count_nonmiss.Consequence.unique()) > 6:
                        ncol = 3
                    else:
                        ncol = 2
                    i = 0
                    axes[ax].vlines(mut_count_nonmiss["Pos"], ymin=0, ymax=mut_count_nonmiss["Count"], color="gray", lw=0.7, zorder=0, alpha=0.5)
                    for cnsq in mut_count_nonmiss.Consequence.unique():
                        count_cnsq = mut_count_nonmiss[mut_count_nonmiss["Consequence"] == cnsq]
                        if cnsq == "synonymous_variant":
                            order = 1
                        else:
                            order = 2
                        if cnsq in plot_pars["color_cnsq"]:
                            color = plot_pars["color_cnsq"][cnsq]
                        else:
                            color=sns.color_palette("tab10")[i]
                            i+=1
                        axes[ax].scatter(count_cnsq.Pos.values, count_cnsq.Count.values, label=capitalize(cnsq), 
                                        color=color, zorder=order, alpha=0.7, lw=plot_pars["s_lw"])              # ec="black",
                    axes[ax].legend(fontsize=11.5, ncol=ncol, framealpha=0.75)
                    axes[ax].set_ylabel('Non\nmissense\nmutations', fontsize=13.5)
                    axes[ax].set_ylim(-0.5, mut_count_nonmiss["Count"].max()+0.5)
                except Exception as e:
                    logger.warning(f"Error occurred while adding non-missense count in {gene} ({uni_id}-F{af_f}):")
                    logger.warning(f"{e}")
                    plot_annot_gene["nonmiss_count"] = False
                    ax -= 1
            else:
                ax -= 1

            # Plot for Missense Mut_in_res track
            # ----------------------------------
            axes[ax+1].vlines(mut_count["Pos"], ymin=0, ymax=mut_count["Count"], color="gray", lw=0.7, zorder=1, alpha=0.5)
            axes[ax+1].scatter(pos_hit, pos_hit_mut, label="Significant", color = 'C0', zorder=4, alpha=0.7, lw=plot_pars["s_lw"])               # ec="black",
            axes[ax+1].scatter(pos_ext, pos_ext_mut, label="Significant extended", color = 'C2', zorder=3, alpha=0.7, lw=plot_pars["s_lw"])      # ec="black",
            axes[ax+1].scatter(pos_not, pos_not_mut, label="Not significant", color = 'C1', zorder=2, alpha=0.7, lw=plot_pars["s_lw"])           # ec="black",
            axes[ax+1].fill_between(pos_result_gene['Pos'], 0, max_mut, where=(pos_result_gene['C'] == 1), 
                            color='skyblue', alpha=0.3, label='Mutated *', zorder=0)
            axes[ax+1].legend(fontsize=11.5, ncol=3, framealpha=0.75)
            axes[ax+1].set_ylabel('Missense\nmutations', fontsize=13.5)         

            # Plot for Score and Miss prob track
            # ----------------------------------
            max_value = np.max((np.max(prob_vec), np.max(score_norm_vec)))
            axes[ax+2].fill_between(pos_result_gene['Pos'], 0, max_value, where=(pos_result_gene['C'] == 1), 
                            color='white')
            axes[ax+2].fill_between(pos_result_gene['Pos'], 0, max_value, where=(pos_result_gene['C'] == 1), 
                            color='skyblue', alpha=0.4, label='Mutated *')
            axes[ax+2].fill_between(range(1, len(prob_vec)+1), 0, score_norm_vec, zorder=1, color="white")
            axes[ax+2].fill_between(range(1, len(prob_vec)+1), 0, score_norm_vec, zorder=1, color="C2", alpha=0.5)            
            axes[ax+2].plot(range(1, len(prob_vec)+1), prob_vec, label="Miss mut prob", zorder=3, color="Red", lw=1)                          
            axes[ax+2].plot(range(1, len(prob_vec)+1), score_norm_vec, label="O3D score normalized", zorder=2, color="C2", lw=0.5)        
            handles, labels = axes[ax+2].get_legend_handles_labels()
            axes[ax+2].legend(handles[-2:], labels[-2:], fontsize=11.5, framealpha=0.75, ncol=2)
            axes[ax+2].set_ylabel('Value', fontsize=13.5)

            # Plot annotations
            # ================

            # Plot PAE
            # ----------------------------------
            if not np.isnan(pos_result_gene["PAE_vol"]).all() and plot_annot_gene["pae"] == True:  
                try:
                    max_value = np.max(pos_result_gene["PAE_vol"])
                    axes[ax+3].fill_between(pos_result_gene['Pos'], 0, max_value, where=(pos_result_gene['C'] == 1), 
                                    color='white')
                    axes[ax+3].fill_between(pos_result_gene['Pos'], 0, max_value, where=(pos_result_gene['C'] == 1), 
                                    color='skyblue', alpha=0.3)
                    axes[ax+3].fill_between(pos_result_gene["Pos"], 0, pos_result_gene["PAE_vol"].fillna(0), 
                                            zorder=2, color=sns.color_palette("pastel")[4])                                           
                    axes[ax+3].plot(pos_result_gene['Pos'], pos_result_gene["PAE_vol"].fillna(0),                                     
                                    label="Confidence", zorder=3, color=sns.color_palette("tab10")[4], lw=0.5)
                    axes[ax+3].set_ylabel('PAE', fontsize=13.5)
                except Exception as e:
                    logger.warning(f"Error occurred while adding PAE in {gene} ({uni_id}-F{af_f}):")
                    logger.warning(f"{e}")
                    plot_aplot_annot_genennot["pae"] = False
                    ax-=1
            else:
                ax-=1
                
            # Plot disorder
            # -------------
            if plot_annot_gene["disorder"]:
                try:
                    axes[ax+4].fill_between(pos_result_gene['Pos'], 0, 100, where=(pos_result_gene['C'] == 1), 
                                            color='white')
                    axes[ax+4].fill_between(pos_result_gene['Pos'], 0, 100, where=(pos_result_gene['C'] == 1), 
                                            color='skyblue', alpha=0.4, label='Mutated *')
                    
                    # ## Comment out to use AF color scale
                    
                    # af_colors = ["#1F6AD7", 
                    #             "#65CBF3",
                    #             "#FFDC48",
                    #             "#FB7C44"]

                    # disorder_x, disorder_y = interpolate_x_y(disorder_gene["Pos"], disorder_gene["Confidence"])
                    # condition_1 = disorder_y > 90
                    # condition_2 = disorder_y <= 90
                    # condition_3 = disorder_y <= 70
                    # condition_4 = disorder_y <= 50
                    # conditions = [condition_1, condition_2, condition_3, condition_4]
                    # for color, condition in zip(af_colors, conditions):
                    #     axes[ax+4].fill_between(disorder_x, 0, disorder_y, where=(condition),       
                    #                             zorder=2, color="white")   
                    #     axes[ax+4].fill_between(disorder_x, 0, disorder_y, where=(condition),   
                    #                             zorder=3, facecolor=color, alpha=0.8)  
                    
                    axes[ax+4].fill_between(disorder_gene["Pos"], 0, disorder_gene["Confidence"].fillna(0),                  
                                zorder=2, color=sns.color_palette("pastel")[4])
                    axes[ax+4].plot(disorder_gene["Pos"], disorder_gene["Confidence"], 
                                    label="Confidence", zorder=3, color=sns.color_palette("tab10")[4], lw=0.5)    
                    axes[ax+4].set_ylabel('pLDDT', fontsize=13.5)
                    axes[ax+4].set_ylim(-10, 110)
                except Exception as e:
                    logger.warning(f"Error occurred while adding Disorder in {gene} ({uni_id}-F{af_f}):")
                    logger.warning(f"{e}")
                    plot_annot_gene["disorder"] = False
                    ax-=1
            else:
                ax-=1

            # Plot pACC
            # ---------
            if plot_annot_gene["pacc"]:
                try:
                    axes[ax+5].fill_between(pos_result_gene['Pos'], 0, 100, where=(pos_result_gene['C'] == 1), 
                                            color='white')
                    axes[ax+5].fill_between(pos_result_gene['Pos'], 0, 100, where=(pos_result_gene['C'] == 1), 
                                            color='skyblue', alpha=0.4)
                    axes[ax+5].fill_between(pdb_tool_gene["Pos"], 0, pdb_tool_gene["pACC"].fillna(0),                  
                                            zorder=2, color=sns.color_palette("pastel")[4])
                    axes[ax+5].plot(pdb_tool_gene['Pos'], pdb_tool_gene["pACC"].fillna(0), 
                                    label="pACC", zorder=3, color=sns.color_palette("tab10")[4], lw=0.5)      
                    axes[ax+5].set_ylabel('pACC', fontsize=13.5)
                    axes[ax+5].set_ylim(-10, 110)
                except Exception as e:
                    logger.warning(f"Error occurred while adding pACC in {gene} ({uni_id}-F{af_f}):")
                    logger.warning(f"{e}")
                    plot_annot_gene["pacc"] = False
            else:
                ax-=1

            # Plot stability change
            # ---------------------
            if pos_result_gene["DDG"].isna().all():
                plot_annot_gene["ddg"] = False
                logger.warning(f"DDG not available for {gene} ({uni_id}-F{af_f}): The track will not be included...")

            if plot_annot_gene["ddg"]:
                try:
                    max_value, min_value = pos_result_gene["DDG"].max(), pos_result_gene["DDG"].min()
                    axes[ax+6].fill_between(pos_result_gene['Pos'], min_value, max_value, where=(pos_result_gene['C'] == 1), 
                                            color='white')
                    axes[ax+6].fill_between(pos_result_gene['Pos'], min_value, max_value, where=(pos_result_gene['C'] == 1), 
                                            color='skyblue', alpha=0.4)
                    axes[ax+6].fill_between(pos_result_gene['Pos'], 0, pos_result_gene["DDG"], zorder=1,             
                                            color=sns.color_palette("pastel")[4])      
                    axes[ax+6].plot(pos_result_gene['Pos'], pos_result_gene["DDG"], 
                                    label="Stability change", zorder=2, color=sns.color_palette("tab10")[4], lw=0.5)    
                    axes[ax+6].set_ylabel('DDG', fontsize=13.5)
                except Exception as e:
                    logger.warning(f"Error occurred while adding DDG in {gene} ({uni_id}-F{af_f}):")
                    logger.warning(f"{e}")
                    plot_annot_gene["ddg"] = False
            else:
                ax-=1 

            # PTM
            # --------------
            if plot_annot_gene["ptm"]:    
                try:
                    ptm_gene = uni_feat_gene[uni_feat_gene["Type"] == "PTM"]
                    ptm_names = ptm_gene["Description"].unique()
                    sb_width = 0.5
                    max_value = (len(ptm_names) * sb_width) - 0.2
                    min_value = - 0.3
                    
                    axes[ax+7].fill_between(pos_result_gene['Pos'], min_value, max_value, where=(pos_result_gene['C'] == 1), 
                                            color='white')
                    axes[ax+7].fill_between(pos_result_gene['Pos'], min_value, max_value, where=(pos_result_gene['C'] == 1), 
                                            color='skyblue', alpha=0.4, label='Mutated *')

                    for n, name in enumerate(ptm_names):
                        c = sns.color_palette("tab10")[n]
                        ptm = ptm_gene[ptm_gene["Description"] == name]
                        ptm_pos = ptm.Begin.values
                        axes[ax+7].scatter(ptm_pos, np.repeat(n*sb_width, len(ptm_pos)), label=name, alpha=0.7, color=c) #label=name
                        axes[ax+7].hlines(y=n*sb_width, xmin=0, xmax=gene_len, linewidth=1, color='lightgray', alpha=0.7, zorder=0)
                
                    axes[ax+7].set_ylim(min_value, max_value)
                    axes[ax+7].set_yticks(sb_width * np.arange(len(ptm_names)), ptm_names)
                    axes[ax+7].set_ylabel(' PTM            ', fontsize=13.5, rotation=0, va='center')
                except Exception as e:
                    logger.warning(f"Error occurred while adding PTM in {gene} ({uni_id}-F{af_f}):")
                    logger.warning(f"{e}")
                    plot_annot_gene["ptm"] = False
            else:
                ax-=1

            # SITES
            # --------------
            if plot_annot_gene["site"]:   
                try:
                    site_gene = uni_feat_gene[uni_feat_gene["Type"] == "SITE"]
                    site_names = site_gene["Description"].unique()
                    sb_width = 0.5
                    max_value = (len(site_names) * sb_width) - 0.2
                    min_value = - 0.3
                    
                    axes[ax+8].fill_between(pos_result_gene['Pos'], min_value, max_value, where=(pos_result_gene['C'] == 1), 
                                            color='white')
                    axes[ax+8].fill_between(pos_result_gene['Pos'], min_value, max_value, where=(pos_result_gene['C'] == 1), 
                                            color='skyblue', alpha=0.4, label='Mutated *')
                    
                    for n, name in enumerate(site_names):
                        c = sns.color_palette("tab10")[n]
                        site = site_gene[site_gene["Description"] == name]
                        site_pos = site.Begin.values
                        axes[ax+8].scatter(site_pos, np.repeat(n*sb_width, len(site_pos)), label=name, alpha=0.7, color=c) #label=name
                        axes[ax+8].hlines(y=n*sb_width, xmin=0, xmax=gene_len, linewidth=1, color='lightgray', alpha=0.7, zorder=0)
                    
                    axes[ax+8].set_ylim(min_value, max_value)
                    axes[ax+8].set_yticks(sb_width * np.arange(len(site_names)), site_names)
                    axes[ax+8].set_ylabel('Site           ', fontsize=13.5, rotation=0, va='center')
                except Exception as e:
                    logger.warning(f"Error occurred while adding Sites in {gene} ({uni_id}-F{af_f}):")
                    logger.warning(f"{e}")
                    plot_annot_gene["site"] = False
            else:
                ax-=1

            # Clusters label
            # --------------
            if plot_annot_gene["clusters"]:
                try:
                    clusters_label = pos_result_gene.Cluster.dropna().unique()
                    palette = sns.color_palette(cc.glasbey, n_colors=len(clusters_label))
                    for i, cluster in enumerate(clusters_label):
                        axes[ax+9].fill_between(pos_result_gene['Pos'], -0.5, 0.46, 
                                                where=((pos_result_gene['Cluster'] == cluster) & (pos_result_gene['C'] == 1)),
                                                color=palette[i], lw=0.4) # alpha=0.6
                    axes[ax+9].set_ylabel('Clusters             ', fontsize=13.5, rotation=0, va='center')
                    axes[ax+9].set_yticks([])  
                    axes[ax+9
                    ].set_yticklabels([], fontsize=12)
                except Exception as e:
                    logger.warning(f"Error occurred while adding Clusters labels in {gene} ({uni_id}-F{af_f}):")
                    logger.warning(f"{e}")
                    plot_annot_gene["clusters"] = False
            else:
                ax-=1
                
            # Secondary structure
            # -------------------
            if plot_annot_gene["sse"]:
                try:
                    for i, sse in enumerate(['Helix', 'Ladder', 'Coil']):
                        c = 0+i
                        ya, yb = c-plot_pars["sse_fill_width"], c+plot_pars["sse_fill_width"]
                        axes[ax+10].fill_between(pdb_tool_gene["Pos"].values, ya, yb, where=(pdb_tool_gene["SSE"] == sse), 
                                        color=sns.color_palette("tab10")[7+i], label=sse)
                    axes[ax+10].set_yticks([0, 1, 2])  
                    axes[ax+10].set_yticklabels(['Helix', 'Ladder', 'Coil'], fontsize=10)
                    axes[ax+10].set_ylabel(' SSE', fontsize=13.5)
                except Exception as e:
                    logger.warning(f"Error occurred while adding SSE in {gene} ({uni_id}-F{af_f}):")
                    logger.warning(f"{e}")
                    plot_annot_gene["sse"] = False
            else:
                ax-=1

            # Pfam
            # ----
            if plot_annot_gene["pfam"]:
                pfam_gene = uni_feat_gene[(uni_feat_gene["Type"] == "DOMAIN") & (uni_feat_gene["Evidence"] == "Pfam")]
                try:
                    pfam_gene = pfam_gene.sort_values("Begin").reset_index(drop=True)
                    pfam_color_dict = {}
                    
                    for n, name in enumerate(pfam_gene["Description"].unique()):
                        pfam_color_dict[name] = f"C{n}"
                        
                    n = 0
                    added_pfam = []
                    for i, row in pfam_gene.iterrows():
                        if pd.Series([row["Description"], row["Begin"], row["End"]]).isnull().any():
                            continue
                        
                        name = row["Description"]
                        start = int(row["Begin"])
                        end = int(row["End"])
                        axes[ax+11].fill_between(range(start, end+1), -0.45, 0.45,  alpha=0.5, color=pfam_color_dict[name])
                        if name not in added_pfam:
                            if near_pfam:
                                n += 1
                                if n == 1:
                                    y = 0.28
                                elif n == 2:
                                    y = 0
                                elif n == 3:
                                    y = -0.295
                                    n = 0
                            else:
                                y = -0.04
                            axes[ax+11].text(((start + end) / 2)+0.5, y, name, ha='center', va='center', fontsize=10, color="black")
                            added_pfam.append(name)
                    axes[ax+11].set_yticks([])  
                    axes[ax+11].set_yticklabels([], fontsize=12)
                    axes[ax+11].set_ylabel('Pfam        ', fontsize=13.5, rotation=0, va='center')
                    axes[ax+11].set_ylim(-0.5, 0.5)  
                except Exception as e:
                    logger.warning(f"Error occurred while adding Pfam domain in {gene} ({uni_id}-F{af_f}):")
                    logger.warning(f"{e}")
                    plot_annot_gene["pfam"] = False
            else:
                ax-=1

            # Prosite
            # -------
            if plot_annot_gene["prosite"]:
                try:
                    prosite_gene = uni_feat_gene[(uni_feat_gene["Type"] == "DOMAIN") & (uni_feat_gene["Evidence"] != "Pfam")]
                
                    prosite_gene = prosite_gene.sort_values("Begin").reset_index(drop=True)
                    prosite_color_dict = {}
                    
                    for n, name in enumerate(prosite_gene["Description"].unique()):
                        prosite_color_dict[name] = f"C{n}"
                        
                    n = 0
                    added_prosite = []
                    for i, row in prosite_gene.iterrows():
                        if pd.Series([row["Description"], row["Begin"], row["End"]]).isnull().any():
                            continue
                        
                        name = row["Description"]
                        start = int(row["Begin"])
                        end = int(row["End"])
                        axes[ax+12].fill_between(range(start, end+1), -0.45, 0.45,  alpha=0.5, color=prosite_color_dict[name])
                        if name not in added_prosite:
                            if near_prosite:
                                n += 1
                                if n == 1:
                                    y = 0.28
                                elif n == 2:
                                    y = 0
                                elif n == 3:
                                    y = -0.295
                                    n = 0
                            else:
                                y = -0.04
                            axes[ax+12].text(((start + end) / 2)+0.5, y, name, ha='center', va='center', fontsize=10, color="black")
                            added_prosite.append(name)
                    axes[ax+12].set_yticks([])  
                    axes[ax+12].set_yticklabels([], fontsize=12)
                    axes[ax+12].set_ylabel('Prosite           ', fontsize=13.5, rotation=0, va='center')
                    axes[ax+12].set_ylim(-0.5, 0.5)  
                except Exception as e:
                    logger.warning(f"Error occurred while adding Prosite domain in {gene} ({uni_id}-F{af_f}):")
                    logger.warning(f"{e}")
                    plot_annot_gene["pfam"] = False
            else:
                ax-=1

            # Membrane
            # --------
            if plot_annot_gene["membrane"]:
                try:
                    membrane_gene = uni_feat_gene[(uni_feat_gene["Type"] == "MEMBRANE")]
                    
                    membrane_gene = membrane_gene.sort_values("Begin").reset_index(drop=True)
                    membrane_color_dict = {}
                    
                    for n, name in enumerate(membrane_gene["Description"].unique()):
                        membrane_color_dict[name] = f"C{n}"
                        
                    n = 0
                    added_membrane = []
                    for i, row in membrane_gene.iterrows():
                        if pd.Series([row["Description"], row["Begin"], row["End"]]).isnull().any():
                            continue
                        
                        name = row["Description"]
                        start = int(row["Begin"])
                        end = int(row["End"])
                        axes[ax+14].fill_between(range(start, end+1), -0.45, 0.45,  alpha=0.5, color=membrane_color_dict[name])
                        if name not in added_membrane:
                            y = -0.04
                            axes[ax+14].text(((start + end) / 2)+0.5, y, name, ha='center', va='center', fontsize=10, color="black")
                            added_membrane.append(name)
                    axes[ax+13].set_yticks([])  
                    axes[ax+13].set_yticklabels([], fontsize=12)
                    axes[ax+13].set_ylabel('Membrane                ', fontsize=13.5, rotation=0, va='center')
                    axes[ax+13].set_ylim(-0.5, 0.5)  
                except Exception as e:
                    logger.warning(f"Error occurred while adding Membrane info in {gene} ({uni_id}-F{af_f}):")
                    logger.warning(f"{e}")
                    plot_annot_gene["membrane"] = False
            else:
                ax-=1

            # Motifs
            # ------
            if plot_annot_gene["motif"]:
                try:
                    motif_gene = uni_feat_gene[(uni_feat_gene["Type"] == "MOTIF")]
                    
                    motif_gene = motif_gene.sort_values("Begin").reset_index(drop=True)
                    motif_color_dict = {}
                    
                    for n, name in enumerate(motif_gene["Full_description"].unique()):
                        motif_color_dict[name] = f"C{n}"
                        
                    n = 0
                    added_motif = []
                    for i, row in motif_gene.iterrows():
                        if pd.Series([row["Full_description"], row["Begin"], row["End"]]).isnull().any():
                            continue
                        
                        name = row["Full_description"]
                        start = int(row["Begin"])
                        end = int(row["End"])
                        axes[ax+14].fill_between(range(start, end+1), -0.45, 0.45,  alpha=0.5, color=motif_color_dict[name])
                        if name not in added_motif:
                            if near_motif:
                                n += 1
                                if n == 1:
                                    y = 0.28
                                elif n == 2:
                                    y = 0
                                elif n == 3:
                                    y = -0.295
                                    n = 0
                            else:
                                y = -0.04
                            axes[ax+14].text(((start + end) / 2)+0.5, y, name, ha='center', va='center', fontsize=10, color="black")
                            added_motif.append(name)
                    axes[ax+14].set_yticks([])  
                    axes[ax+14].set_yticklabels([], fontsize=12)
                    axes[ax+14].set_ylabel('Motif        ', fontsize=13.5, rotation=0, va='center')
                    axes[ax+14].set_ylim(-0.5, 0.5) 
                except Exception as e:
                    logger.warning(f"Error occurred while adding Motifs info in {gene} ({uni_id}-F{af_f}):")
                    logger.warning(f"{e}")
                    plot_annot_gene["motif"] = False
            else:
                ax-=1

            # Save
            # ----
            fig.suptitle(f'{gene} - {uni_id}', fontsize=16)
            filename = f"{cohort}.genes_plot_{j+1}.{gene}_{uni_id}.png"
            output_path = os.path.join(output_dir, filename)
            plt.subplots_adjust(top=0.95) 
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            logger.debug(f"Saved {output_path}")
            plt.close()
            
            # Store annotated result
            pos_result_gene = get_enriched_result(pos_result_gene, 
                                                  disorder_gene, 
                                                  pdb_tool_gene, 
                                                  seq_df)
            annotated_result_lst.append(pos_result_gene)
            uni_feat_result_lst.append(uni_feat_gene)

    # Save as tsv
    if output_tsv:
        pos_result_annotated = pd.concat(annotated_result_lst)
        feat_processed = pd.concat(uni_feat_result_lst)   
    else:
        pos_result_annotated = None
        feat_processed = None
        
    return pos_result_annotated, feat_processed
        
            
# ============
# PLOT WRAPPER
# ============

def generate_plot(gene_result_path, 
                  pos_result_path, 
                  maf_path, 
                  mut_profile_path,
                  mutability_config_path,
                  gene_result_path_2,
                  pos_result_path_2,
                  maf_path_2,
                  mut_profile_path_2,
                  mutability_config_path_2,
                  datasets_dir, 
                  annotations_dir,
                  output_dir,
                  cohort,
                  plot_annot,
                  plot_pars,
                  n_genes=30, 
                  non_significant=False, 
                  lst_genes=None,
                  comparative_plots=False,
                  output_tsv=False,
                  output_all_pos=False):
    
    # Load data tracks
    # ================
    
    # Load data
    gene_result = pd.read_csv(gene_result_path, sep="\t")
    pos_result = pd.read_csv(pos_result_path, sep="\t")
    maf = parse_maf_input(maf_path)
    
    seq_df_path = os.path.join(datasets_dir, "seq_for_mut_prob.tsv") 
    seq_df = pd.read_csv(seq_df_path, sep="\t")    
    uniprot_feat = pd.read_csv(os.path.join(annotations_dir, "uniprot_feat.tsv"), sep="\t")
    uniprot_feat.loc[(uniprot_feat["Type"] == "MOTIF") & (
        uniprot_feat["Description"] == "Zinc finger"), "Full_description"] = "Zinc finger"
    pdb_tool = pd.read_csv(os.path.join(annotations_dir, "pdb_tool_df.tsv"), sep="\t")
    disorder = pd.read_csv(os.path.join(datasets_dir, "confidence.tsv"), sep="\t", low_memory=False)
    dict_transcripts = plot_pars["dict_transcripts"]

    # Filter Oncodrive3D result
    gene_result, pos_result, genes, uni_ids = filter_o3d_result(gene_result, 
                                                                pos_result, 
                                                                n_genes, 
                                                                lst_genes, 
                                                                non_significant)
    
    if len(gene_result) > 0:   

        # Subset dfs by selected genes and IDs
        seq_df, disorder, pdb_tool, uniprot_feat = subset_genes_and_ids(genes, 
                                                                        uni_ids, 
                                                                        seq_df, 
                                                                        disorder, 
                                                                        pdb_tool, 
                                                                        uniprot_feat)

        # Get missense mut probability dict
        miss_prob_dict = get_miss_mut_prob_for_plot(mut_profile_path, mutability_config_path, seq_df)

        # Summary plots
        # =============
        
        create_plot_dir(output_dir)
        logger.info(f"Creating summary plot in {output_dir}")
        count_mut_gene_df, count_pos_df, cluster_df, size_df = get_summary_counts(gene_result, pos_result, seq_df)
        summary_plot(gene_result, 
                     pos_result, 
                     count_mut_gene_df, 
                     count_pos_df, 
                     cluster_df,
                     size_df,
                     output_dir,
                     cohort) 
        
        # Plots for individual genes
        # ==========================
        
        # Get non missense mutations
        if plot_annot["nonmiss_count"]:
            maf_nonmiss = get_nonmiss_mut(maf_path)
        else:
            maf_nonmiss = None
        if maf_nonmiss is None:
            plot_annot["nonmiss_count"] = False
        
        output_dir_genes_plots = os.path.join(output_dir, f"{cohort}.genes_plots")
        create_plot_dir(output_dir_genes_plots)
        logger.info(f"Creating genes plots in {output_dir_genes_plots}")
        pos_result_annotated, uni_feat_processed = genes_plots(gene_result, 
                                                            pos_result, 
                                                            seq_df,
                                                            maf,
                                                            miss_prob_dict,
                                                            output_dir_genes_plots,
                                                            cohort,
                                                            annotations_dir,
                                                            disorder,
                                                            uniprot_feat,
                                                            pdb_tool,
                                                            maf_nonmiss,
                                                            plot_annot,
                                                            plot_pars,
                                                            size_df,
                                                            output_tsv)
        
        # Save annotations
        if output_tsv and pos_result_annotated is not None:
            logger.info(f"Saving annotated Oncodrive3D result in {output_dir}")
            save_annotated_pos_result(pos_result, 
                                       pos_result_annotated, 
                                       uni_feat_processed, 
                                       output_dir, 
                                       cohort, 
                                       output_all_pos)
            
        logger.info("Plotting completed!")
    
    else:
        logger.warning("There aren't any genes to plot!")
        
        
if __name__ == "__main__":  
    
    from scripts.plotting.utils import init_annotations
    
    # Args
    cohort = "PEDCBIOP_WXS_OS_PRY"
    #cohort = "CPTAC_WXS_BRCA_2020"
    figsize_x=24
    figsize_y=12
    
    # Annotations
    plot_annot = init_annotations("all")

    # Cnsq color
    color_cnsq = {"splicing" : "C2",
                  "missense" : "C5",
                  "synonymous" : "C9",
                  "coding_sequence_variant" : "C1",
                  "nonsense" : "C6",
                  "intron_variant" : "C7",
                  "indel" : "C8",
                  "protein_altering_variant" : "C3"}

    # Plot parameters                                            # TODO: add some of them as args
    plot_pars = {}
    plot_pars["figsize"] = figsize_x, figsize_y
    plot_pars["h_ratios"] = [0.17, 0.17, 0.17, 0.1, 0.1, 0.1, 0.1, 0.022, 0.022, 0.04, 0.07, 0.04, 0.04, 0.04, 0.04]
    plot_pars["s_lw"] = 0.2
    plot_pars["sse_fill_width"] = 0.43
    plot_pars["dist_thr"] = 0.05
    plot_pars["color_cnsq"] = color_cnsq
    plot_pars["dict_transcripts"] = {"PTEN" : "ENST00000688308"}  # For now this is used to get a specific transcript for pfam domains
    
    path_o3d_result = "/workspace/projects/clustering_3d/o3d_analysys/datasets/output/cancer/o3d_output/run_2024-03-06_11-47-05"
    path_o3d_result = "/workspace/projects/clustering_3d/o3d_analysys/datasets/output/cancer/o3d_output/run_2024-03-19_10-52-36"
    path_input = "/workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer/"
    
    generate_plot(gene_result_path=f"{path_o3d_result}/{cohort}/{cohort}.3d_clustering_genes.tsv",
                  pos_result_path=f"{path_o3d_result}/{cohort}/{cohort}.3d_clustering_pos.tsv",
                  maf_path=f"{path_input}/maf/{cohort}.in.maf",
                  mut_profile_path=f"{path_input}/mut_profile/{cohort}.mutrate.json",
                  mutability_config_path=None,
                  gene_result_path_2=None,
                  pos_result_path_2=None,
                  maf_path_2=None,
                  mut_profile_path_2=None,
                  mutability_config_path_2=None,
                  datasets_dir="/workspace/nobackup/scratch/oncodrive3d/datasets_mane",
                  annotations_dir="/workspace/nobackup/scratch/oncodrive3d/annotations_mane",
                  output_dir=f"{path_o3d_result}/{cohort}",
                  cohort=cohort,
                  plot_annot=plot_annot,
                  plot_pars=plot_pars,
                  n_genes=30,
                  lst_genes=None,
                  non_significant=True,
                  comparative_plots=False,
                  output_tsv=True,
                  output_all_pos=False)