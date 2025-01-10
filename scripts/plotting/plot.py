import pandas as pd
import daiquiri
import logging
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import os
import json
import colorcet as cc
import warnings
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import ConvergenceWarning
from adjustText import adjust_text
from matplotlib.axes._axes import Axes
from scripts.plotting.utils import get_broad_consequence, save_annotated_result
from scripts.plotting.utils import get_enriched_result, filter_o3d_result, subset_genes_and_ids, load_o3d_result

from scripts import __logger_name__

logger = daiquiri.getLogger(__logger_name__ + ".plotting.plot")

logging.getLogger('matplotlib.font_manager').setLevel(logging.WARNING)



# Summary plots
# =============

def get_summary_counts(gene_result, pos_result, seq_df):
    """
    Get dataframes including the counts required to generate the summary plots.
    """

    # Df with mut count
    count_mut_gene_hit = gene_result[["Gene", "Clust_mut"]].rename(columns={"Clust_mut" : "Mut_in_gene"})
    count_mut_gene_hit["C"] = "Mutations in clusters"
    count_mut_gene_not = pd.DataFrame({"Gene" : count_mut_gene_hit["Gene"].values,
                                    "Mut_in_gene" : gene_result.apply(lambda x: x["Mut_in_gene"] - x["Clust_mut"], axis=1)})
    count_mut_gene_not["C"] = "Mutations not in clusters"
    count_mut_gene = gene_result[["Gene", "Mut_in_gene"]].copy()
    count_mut_gene["C"] = "Total mutations"

    count_mut_gene_df = pd.concat((count_mut_gene_hit, count_mut_gene_not, count_mut_gene)).sort_values("Gene").rename(columns={"Mut_in_gene" : "Count"})
    count_mut_gene_df = count_mut_gene_df.sort_values(["C", "Count"], ascending=False).reset_index(drop=True)

    # Df with pos count
    pos_result_not = pos_result[pos_result["C"] == 0]
    if len(pos_result_not) > 0: 
        pos_result_not = pos_result_not.groupby("Gene").apply(len)
        pos_result_not = pos_result_not.reset_index().rename(columns={0 : "Count"})
        pos_result_not["C"] = "Residues not in clusters"
    else:
        pos_result_not = pd.DataFrame(columns=["Gene", "Count", "C"])
    pos_result_hit = pos_result[pos_result["C"] == 1]
    if len(pos_result_hit) > 0:   
        pos_result_hit = pos_result_hit.groupby("Gene").apply(len)
        pos_result_hit = pos_result_hit.reset_index().rename(columns={0 : "Count"})
        pos_result_hit["C"] = "Residues in clusters"
    else:
        pos_result_hit = pd.DataFrame(columns=["Gene", "Count", "C"])

    pos_result_total = pd.DataFrame(seq_df.apply(lambda x: (x.Gene, len(x.Seq)), axis=1).to_list())
    pos_result_total.columns = "Gene", "Count"
    pos_result_total["C"] = "Protein length"
    
    count_pos_df = pd.concat((pos_result_total, pos_result_hit, pos_result_not)).sort_values("Gene")
    count_pos_df = count_pos_df.sort_values("C", ascending=False).reset_index(drop=True)

    # Df with cluster count
    cluster_df = pos_result.groupby("Gene").max("Clump").Clump.reset_index()
    cluster_df["Clump"] = cluster_df["Clump"] + 1
    
    return count_mut_gene_df, count_pos_df, cluster_df
        

def summary_plot(gene_result, 
                 pos_result, 
                 count_mut_gene_df, 
                 count_pos_df, 
                 cluster_df,
                 output_dir,
                 cohort,
                 plot_pars,
                 save_plot=True,
                 show_plot=False,
                 title=None):        

    # Init
    h_ratios = plot_pars["summary_h_ratios"]
    tracks = list(h_ratios.keys())
    
    # Plot
    fsize_x, fsize_y = plot_pars["summary_figsize"]
    if len(gene_result) < 6:
        fsize_x = 3
    else:
        fsize_x = fsize_x * len(gene_result)
    fig, axes = plt.subplots(len(h_ratios), 1, 
                             figsize=(fsize_x, fsize_y), 
                             sharex=True, 
                             gridspec_kw={'hspace': 0.1, 
                                          'height_ratios': h_ratios.values()})
    
    if "score" in tracks:
        ax = tracks.index("score")
        pos_result = pos_result.copy()
        pos_result["C"] = pos_result.C.map({1 : "Volume in clusters", 0 : "Volume not in clusters"})
        hue_order = ['Volume in clusters', 'Volume not in clusters']
        sns.boxplot(x='Gene', y='Score_obs_sim', data=pos_result, order=gene_result.Gene, color=sns.color_palette("pastel")[7], showfliers=False, ax=axes[ax])
        sns.stripplot(x='Gene', y='Score_obs_sim', data=pos_result, hue="C" ,jitter=True, size=6, alpha=plot_pars["summary_alpha"], order=gene_result.Gene.values, 
                      palette=sns.color_palette("tab10", n_colors=2), hue_order=hue_order, ax=axes[ax])
        axes[ax].set_ylabel('Clustering\nscore\n(obs/sim)', fontsize=12)
        axes[ax].legend(fontsize=9.5, loc="upper right")
        axes[ax].set_xlabel(None)
    
    if "miss_count" in tracks:
        ax = tracks.index("miss_count")
        hue_order = ['Total mutations', 'Mutations in clusters']
        custom_palette = [sns.color_palette("pastel")[7], sns.color_palette("pastel")[0]]
        sns.barplot(x='Gene', y='Count', data=count_mut_gene_df[count_mut_gene_df["C"] != "Mutations not in clusters"], 
                    order=gene_result.Gene, ax=axes[ax], hue="C", palette=custom_palette, hue_order=hue_order, ec="black", lw=0.5)
        axes[ax].set_ylabel('Missense\nmut count', fontsize=12)
        axes[ax].legend(fontsize=9.5, loc="upper right")
        axes[ax].set_xlabel(None)
    
    if "res_count" in tracks:
        ax = tracks.index("res_count")
        hue_order = ['Protein length', 'Residues in clusters']
        sns.barplot(x='Gene', y='Count', data=count_pos_df[count_pos_df["C"] != "Residues not in clusters"], order=gene_result.Gene, hue="C", ax=axes[ax],
                    palette=custom_palette, hue_order=hue_order, ec="black", lw=0.5)
        axes[ax].set_ylabel('Residues\ncount', fontsize=12)
        axes[ax].legend(fontsize=9.5, loc="upper right")
        axes[ax].set_xlabel(None)
    
    if "res_clust_mut" in tracks:
        ax = tracks.index("res_clust_mut")
        count_mut_clusters = count_mut_gene_df[count_mut_gene_df["C"] == "Mutations in clusters"].drop(columns="C").rename(columns={"Count" : "Mut_count"})
        count_res_clusters = count_pos_df[count_pos_df["C"] == "Residues in clusters"].drop(columns="C").rename(columns={"Count" : "Res_count"})
        count_mut_res_clusters = count_mut_clusters.merge(count_res_clusters, on="Gene")
        count_mut_res_clusters["Per_res_mut"] = np.round(count_mut_res_clusters["Mut_count"] / count_mut_res_clusters["Res_count"], 2)
        sns.barplot(x='Gene', y='Per_res_mut', data=count_mut_res_clusters, order=gene_result.Gene, ax=axes[ax], color=sns.color_palette("pastel")[0], ec="black", lw=0.5)
        axes[ax].set_ylabel('Per-residue\nmut\nin clusters', fontsize=12)
        axes[ax].set_xlabel(None)
    
    if "clusters" in tracks:
        ax = tracks.index("clusters")
        sns.barplot(x='Gene', y='Clump', data=cluster_df, order=gene_result.Gene, ax=axes[ax], color=sns.color_palette("pastel")[0], ec="black", lw=0.5)
        axes[ax].set_ylabel('Clumps', fontsize=12)
        axes[ax].set_xlabel(None)
    
    # Details
    if title: 
        fig.suptitle(f"{title} summary", fontsize=14)
    else:
        fig.suptitle(f"O3D analysis summary", fontsize=14)
    xticks_labels = [ r'$\mathbf{*}$ ' + gene if gene_result.loc[gene_result["Gene"] == gene, "C_gene"].values[0] == 1 else gene for gene in gene_result.Gene]
    axes[len(axes)-1].set_xticklabels(xticks_labels, rotation=45, rotation_mode="anchor", ha='right', fontsize=12)
    plt.xticks(rotation=45, rotation_mode="anchor", ha='right', fontsize=12)
    plt.subplots_adjust(top=0.94) 
    
    # Save
    filename = f"{cohort}.summary_plot.png"
    output_path = os.path.join(output_dir, filename)
    if save_plot: 
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.debug(f"Saved {output_path}")
    if show_plot:
        plt.show()
    plt.close()


# Gene plots
# ==========

def check_near_feat(uni_feat_gene, feat, dist_thr=0.05):
    """
    Check if two domains could be closer to each other 
    than allowed threshold (ratio of protein size).
    """

    near_feat = False
    uni_feat_gene = uni_feat_gene.copy()
    
    if feat == "domain" or feat == "pfam":
        uni_feat_gene = uni_feat_gene[uni_feat_gene["Type"] == "DOMAIN"]
        if feat == "pfam":
            uni_feat_gene = uni_feat_gene[uni_feat_gene["Evidence"] == "Pfam"]
        else:
            uni_feat_gene = uni_feat_gene[uni_feat_gene["Evidence"] != "Pfam"]
        uni_feat_gene = uni_feat_gene.drop_duplicates(subset='Description', keep='first')
    elif feat == "motif":
        uni_feat_gene = uni_feat_gene[uni_feat_gene["Type"] == "MOTIF"]
        uni_feat_gene = uni_feat_gene.drop_duplicates(subset='Full_description', keep='first')

    mid_pos = (uni_feat_gene.Begin + uni_feat_gene.End) / 2
    mid_pos_norm = (mid_pos / mid_pos.max()).values
    
    for i in range(len(mid_pos_norm)):
        for j in range(i + 1, len(mid_pos_norm)):
            diff = abs(mid_pos_norm[i] - mid_pos_norm[j])
            if diff < dist_thr:
                near_feat = True

    return near_feat


def get_gene_arg(pos_result_gene, plot_pars, uni_feat_gene, maf_nonmiss=None):
    """
    Adjust the height ratio of tracks to include in the plot. 
    """

    h_ratios = plot_pars["h_ratios"].copy()
    plot_pars = plot_pars.copy()
    
    track = "maf_nonmiss"
    if track in h_ratios and not maf_nonmiss:
        del h_ratios[track]
        
    track = "maf_nonmiss_2"
    if track in h_ratios and not maf_nonmiss:
        del h_ratios[track]
        
    track = "ddg"
    if track in h_ratios and pos_result_gene["DDG"].isna().all():
        del h_ratios[track]
    
    track = "pae"
    if track in h_ratios and np.isnan(pos_result_gene["PAE_vol"]).all():
        del h_ratios[track]
        
    track = "ptm" 
    if track in h_ratios:
        if len(uni_feat_gene[uni_feat_gene["Type"] == "PTM"]) == 0:
            del h_ratios[track]
        else:
            stracks = len(uni_feat_gene[uni_feat_gene["Type"] == "PTM"].Description.unique())
            h_ratios[track] = h_ratios[track] * stracks
 
    track = "site"
    if track in h_ratios:
        if len(uni_feat_gene[uni_feat_gene["Type"] == "SITE"]) == 0:
            del h_ratios[track]
        else:
            stracks = len(uni_feat_gene[uni_feat_gene["Type"] == "SITE"].Description.unique())
            h_ratios[track] = h_ratios[track] * stracks

    track = "pfam"
    if track in h_ratios:
        if len(uni_feat_gene[(uni_feat_gene["Type"] == "DOMAIN") & (uni_feat_gene["Evidence"] == "pfam")]) == 1:
            del h_ratios[track]
            near_pfam = False
        else:
            near_pfam = check_near_feat(uni_feat_gene, feat=track, dist_thr=plot_pars["dist_thr"])
            if near_pfam:
                h_ratios[track] = h_ratios[track] * 2
          
    track = "prosite"
    if track in h_ratios:
        if len(uni_feat_gene[(uni_feat_gene["Type"] == "DOMAIN") & (uni_feat_gene["Evidence"] != "Pfam")]) == 0:
            del h_ratios[track]
            near_prosite = False
        else:
            near_prosite = check_near_feat(uni_feat_gene, feat="domain", dist_thr=plot_pars["dist_thr"])
            if near_prosite:
                 h_ratios[track] = h_ratios[track] * 2
        
    track = "membrane"    
    if track in h_ratios:
        if len(uni_feat_gene[uni_feat_gene["Type"] == "MEMBRANE"]) == 0:
            del h_ratios[track]
        else:
            stracks = len(uni_feat_gene[uni_feat_gene["Type"] == "MEMBRANE"].Description.unique())
            h_ratios[track] = h_ratios[track] * stracks
    
    track = "motif"
    if track in h_ratios:
        if len(uni_feat_gene[uni_feat_gene["Type"] == "MOTIF"]) == 0:
            del h_ratios[track]
            near_motif = False
        else:
            near_motif = check_near_feat(uni_feat_gene, feat="motif", dist_thr=0.1)
            if near_motif:
                h_ratios[track] = h_ratios[track] * 1.8

    h_ratios = {k:v/sum(h_ratios.values()) for k,v in h_ratios.items()}
    
    return h_ratios, near_pfam, near_prosite, near_motif


def filter_non_processed_mut(maf, pos_result):
    """
    Get rid of mutations of the input file that were not processed.
    """
                                                # TO DO: In this way, I am not getting rid of the mismatches ones.. I should think about something else for these ones
    len_maf = len(maf)
    maf = maf[maf.apply(lambda x: f"{x.Gene}_{x.Pos}", axis=1).isin(pos_result.apply(lambda x: f"{x.Gene}_{x.Pos}", axis=1))]
    logger.debug(f"Filtered out {len_maf - len(maf)} ({(len_maf - len(maf))/len_maf*100:.2f}%) mutations out of {len_maf} not processed during 3D-clustering analysis!")

    return maf


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
        logger.warning("Can't parse non-missense mutation from MAF file: The track will not be included...")
        logger.warning(f"{e}")


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


def parse_pos_result_for_genes_plot(pos_result_gene, c_ext=True):
    """
    Get mut count and score divided by Oncodriv3D 
    result: significant, not significant, significant extended
    (mutation in a non-significant residue that contribute to 
    mutations in the volume of a significant one/s).
    """
    
    pos_result_gene = pos_result_gene.copy()
    pos_result_gene = pos_result_gene[["Pos", "Mut_in_res", "Mut_in_vol", "Score_obs_sim", "C", "C_ext", "pval", "Clump", "PAE_vol"]]
    if not c_ext:  
        pos_result_gene["C"] = pos_result_gene.apply(
            lambda x: 1 if (x["C"] == 1) & (x["C_ext"] == 0) else 2 if (x["C"] == 1) & (x["C_ext"] == 1) else 0, axis=1)
    max_mut = np.max(pos_result_gene["Mut_in_res"].values)
    
    return pos_result_gene, max_mut
    
    
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
            score = pos_result_gene.loc[pos_result_gene["Pos"] == pos, "Score_obs_sim"].values[0]
        else:
            score = 0
            row_gene = pd.DataFrame({'Pos': [pos], 'Mut_in_res': [0], 'Score_obs_sim': [np.nan], 'C': [np.nan]})
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
        logger.debug(f"Stability change of {uni_id} not found. Path {ddg_path} doesn't exist: Skipping..")

    # Avoid duplicates (Uniprot IDs mapping to the different gene names)
    uni_feat_gene = uni_feat_gene.drop(columns=["Gene", "Ens_Transcr_ID", "Ens_Gene_ID"]).drop_duplicates()
    
    return pos_result_gene, disorder_gene, pdb_tool_gene, uni_feat_gene


def get_site_pos(site_df):

    positions = []
    for begin, end in zip(site_df['Begin'], site_df['End']):
        positions.extend(np.arange(begin, end + 1))
    
    return np.array(positions)


def genes_plots(gene_result, 
                pos_result, 
                seq_df,
                maf,
                maf_nonmiss,
                miss_prob_dict,
                output_dir,
                cohort,
                annotations_dir,
                disorder,
                uniprot_feat,
                pdb_tool,
                plot_pars,
                save_plot=True,
                show_plot=False,
                title=None,
                c_ext=True):   
    """
    Generate a diagnostic plot for each gene showing Oncodrive3D 
    results and annotated features.
    """
    
    annotated_result_lst = []
    uni_feat_result_lst = []
    for j, gene in enumerate(gene_result["Gene"].values):
      
      
        # Load and parse
        # ==============
        
        # IDs
        uni_id = seq_df[seq_df["Gene"] == gene].Uniprot_ID.values[0]
        af_f = seq_df[seq_df["Gene"] == gene].F.values[0]
        gene_len = len(seq_df[seq_df["Gene"] == gene].Seq.values[0])
        maf_gene = maf[maf["Gene"] == gene]
    
        # Parse
        pos_result_gene = pos_result[pos_result["Gene"] == gene].sort_values("Pos").reset_index(drop=True)

        if len(pos_result_gene) > 0:

            pos_result_gene = pos_result_gene[["Pos", "Mut_in_res", "Mut_in_vol", 
                                               "Score_obs_sim", "C", "C_ext", 
                                               "pval", "Clump", "PAE_vol"]]
            pos_result_gene, max_mut = parse_pos_result_for_genes_plot(pos_result_gene, c_ext=c_ext)
            
            # Counts
            mut_count, mut_count_nonmiss = get_count_for_genes_plot(maf_gene, 
                                                                    maf_nonmiss, 
                                                                    gene, 
                                                                    non_missense_count="nonmiss_count" in plot_pars["h_ratios"])
            
            # Get prob vec
            prob_vec = miss_prob_dict[f"{uni_id}-F{af_f}"]                          # TODO: If none, use uniform        <-------------------------- TODO
            
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
            
            h_ratios, near_pfam, near_prosite, near_motif = get_gene_arg(pos_result_gene, plot_pars, uni_feat_gene, maf_nonmiss=maf_nonmiss)
            annotations = list(h_ratios.keys())
            fig, axes = plt.subplots(len(h_ratios), 1, 
                                     figsize=(24,12), 
                                     sharex=True, 
                                     gridspec_kw={'hspace': 0.1, 
                                                  'height_ratios': h_ratios.values()})
            
            
            # Plot for Non-missense mut track   
            # -------------------------------
            if "nonmiss_count" in annotations:
                ax = annotations.index("nonmiss_count")
                if len(mut_count_nonmiss.Consequence.unique()) > 6:
                    ncol = 3
                else:
                    ncol = 2
                i = 0
                axes[ax].vlines(mut_count_nonmiss["Pos"], ymin=0, ymax=mut_count_nonmiss["Count"], 
                                color="gray", lw=0.7, zorder=0, alpha=0.5) # To cover the overlapping needle top part
                axes[ax].scatter(mut_count_nonmiss["Pos"], mut_count_nonmiss["Count"], color='white', zorder=4, lw=plot_pars["s_lw"]) 
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
                                    color=color, zorder=order, alpha=0.7, lw=plot_pars["s_lw"], ec="black")              # ec="black",
                axes[ax].legend(fontsize=11.5, ncol=ncol, framealpha=0.75)
                axes[ax].set_ylabel('Non\nmissense\nmutations', fontsize=13.5, rotation=0, va='center')
                axes[ax].set_ylim(-0.5, mut_count_nonmiss["Count"].max()+0.5)
                axes[ax].set_ylim(0, max(mut_count_nonmiss["Count"])*1.1)
            
            
            # Plot for Missense Mut_in_res track
            # ----------------------------------
            if "miss_count" in annotations:
                ax = annotations.index("miss_count")
            
                axes[ax].vlines(mut_count["Pos"], ymin=0, ymax=mut_count["Count"], color="gray", lw=0.7, zorder=1, alpha=0.5)
                
                mut_pos = pos_result_gene[pos_result_gene["Mut_in_res"] > 0].Pos.values
                mut_res_pos = pos_result_gene[pos_result_gene["Mut_in_res"] > 0].Mut_in_res.values
                # mut_vol_pos = pos_result_gene[pos_result_gene["Mut_in_res"] > 0].Mut_in_vol.values
                
                axes[ax].scatter(mut_pos, mut_res_pos, color='white', zorder=3, lw=plot_pars["s_lw"], ec="white")                  # To cover the overlapping needle top part
                axes[ax].scatter(mut_pos, mut_res_pos, color='gray', zorder=4, alpha=0.7, lw=plot_pars["s_lw"], ec="black", s=60)    
                
                axes[ax].fill_between(pos_result_gene['Pos'], 0, max_mut, where=(pos_result_gene['C'] == 1), 
                                color='skyblue', alpha=0.3, label='Position in cluster', zorder=0, lw=2)
                # axes[ax].fill_between(pos_result_gene['Pos'], 0, max_mut, where=((pos_result_gene["C"] == 0) | (pos_result_gene["C"] == 2)), 
                #                 color='#ffd8b1', alpha=0.6, label='Mutated not *', zorder=0)
                axes[ax].legend(fontsize=11.5, ncol=2, framealpha=0.75)
                axes[ax].set_ylabel('Missense\nmutations', fontsize=13.5, rotation=0, va='center') 
                axes[ax].yaxis.set_label_coords(-0.06, 0.5)
                axes[ax].set_ylim(0-(max(mut_res_pos)*0.04), max(mut_res_pos)*1.1)
                
                
                legend = axes[ax].legend(fontsize=11.5, ncol=1, framealpha=0.75, bbox_to_anchor=(0.97, 1.5), borderaxespad=0.)
                legend.set_title("Global legend")
                legend.get_title().set_fontsize(12)
            
            # Plot for Miss prob track
            # ----------------------------------
            if "miss_prob" in annotations:
                ax = annotations.index("miss_prob")
                
                max_value = np.max(prob_vec)
                axes[ax].fill_between(pos_result_gene['Pos'], 0, max_value, where=(pos_result_gene['C'] == 1), 
                                color='skyblue', alpha=0.3, label='Position in cluster', zorder=0, lw=2)
                
                # axes[ax].hlines(0, xmin=0, xmax=gene_len, color="gray", lw=0.6, zorder=1)
                axes[ax].plot(range(1, len(prob_vec)+1), prob_vec, zorder=3, color="C2", lw=1)                                        
                axes[ax].set_ylabel('Missense\nmut prob', fontsize=13.5, rotation=0, va='center')
                axes[ax].yaxis.set_label_coords(-0.06, 0.5)
            
            # Plot for Score track
            # ----------------------------------
            if "score" in annotations:
                ax = annotations.index("score")
                
                max_value = np.max(score_vec)
                axes[ax].fill_between(pos_result_gene['Pos'], 0, max_value, where=(pos_result_gene['C'] == 1), 
                                color='skyblue', alpha=0.3, label='Position in cluster', zorder=0)
                
                # axes[ax].hlines(0, xmin=0, xmax=gene_len, color="gray", lw=0.7, zorder=1)
                axes[ax].plot(range(1, len(score_vec)+1), score_vec, zorder=2, color="C2", lw=1)                       
                
                # handles, labels = axes[ax].get_legend_handles_labels()
                # axes[ax].legend(fontsize=11.5, framealpha=0.75, ncol=2)
                axes[ax].set_ylabel('Clustering\nscore\n(obs/sim)', fontsize=13.5, rotation=0, va='center')
                axes[ax].yaxis.set_label_coords(-0.06, 0.5)
                
            
            # Plot annotations
            # ================
            
            # Plot PAE
            # ----------------------------------
            if "pae" in annotations:
                ax = annotations.index("pae")
                
                max_value = np.max(pos_result_gene["PAE_vol"])
                # axes[ax+3].fill_between(pos_result_gene['Pos'], 0, max_value, where=((pos_result_gene["C"] == 0) | (pos_result_gene["C"] == 2)), 
                #                 color='#ffd8b1', alpha=0.6)
                axes[ax].fill_between(pos_result_gene['Pos'], 0, max_value, where=(pos_result_gene['C'] == 1), 
                                color='white', lw=2)
                axes[ax].fill_between(pos_result_gene['Pos'], 0, max_value, where=(pos_result_gene['C'] == 1), 
                                color='skyblue', alpha=0.3, lw=2)
                axes[ax].fill_between(pos_result_gene["Pos"], 0, pos_result_gene["PAE_vol"].fillna(0), 
                                        zorder=2, color="white")    
                axes[ax].fill_between(pos_result_gene["Pos"], 0, pos_result_gene["PAE_vol"].fillna(0), 
                                        zorder=2, color=sns.color_palette("pastel")[4], alpha=0.6)    
                axes[ax].plot(pos_result_gene['Pos'], pos_result_gene["PAE_vol"].fillna(0),                                     
                                label="Confidence", zorder=3, color=sns.color_palette("tab10")[4], lw=0.5)
                axes[ax].set_ylabel('Predicted\naligned error\n(Å)', fontsize=13.5, rotation=0, va='center')
                axes[ax].yaxis.set_label_coords(-0.06, 0.5)
            
            # Plot disorder
            # -------------
            if "disorder" in annotations:
                ax = annotations.index("disorder")
            
                axes[ax].fill_between(pos_result_gene['Pos'], 0, 100, where=(pos_result_gene['C'] == 1), 
                                        color='white', lw=2)
                axes[ax].fill_between(pos_result_gene['Pos'], 0, 100, where=(pos_result_gene['C'] == 1), 
                                        color='skyblue', alpha=0.4, label='Mutated *', lw=2)
            
                axes[ax].fill_between(disorder_gene["Pos"], 0, disorder_gene["Confidence"].fillna(0),                  
                                        zorder=2, color="white")
                axes[ax].fill_between(disorder_gene["Pos"], 0, disorder_gene["Confidence"].fillna(0),                  
                                        zorder=2, color=sns.color_palette("pastel")[4], alpha=0.6)
            
                
                axes[ax].plot(disorder_gene["Pos"], disorder_gene["Confidence"], 
                                label="Confidence", zorder=3, color=sns.color_palette("tab10")[4], lw=0.5)    
                axes[ax].set_ylabel('pLDDT\n(disorder)', fontsize=13.5, rotation=0, va='center')
                axes[ax].yaxis.set_label_coords(-0.06, 0.5)
                axes[ax].set_ylim(-10, 110)
            
            # Plot pACC
            # ---------
            if "pacc" in annotations:
                ax = annotations.index("pacc")
                
                # axes[ax+5].fill_between(pos_result_gene['Pos'], 0, 100, where=((pos_result_gene["C"] == 0) | (pos_result_gene["C"] == 2)), 
                #                         color='#ffd8b1', alpha=0.6)
                axes[ax].fill_between(pos_result_gene['Pos'], 0, 100, where=(pos_result_gene['C'] == 1), 
                                        color='white', lw=2)
                axes[ax].fill_between(pos_result_gene['Pos'], 0, 100, where=(pos_result_gene['C'] == 1), 
                                        color='skyblue', alpha=0.4, lw=2)
                axes[ax].fill_between(pdb_tool_gene["Pos"], 0, pdb_tool_gene["pACC"].fillna(0),                  
                                        zorder=2, color="white")
                axes[ax].fill_between(pdb_tool_gene["Pos"], 0, pdb_tool_gene["pACC"].fillna(0),                  
                                        zorder=2, color=sns.color_palette("pastel")[4], alpha=0.6)
                axes[ax].plot(pdb_tool_gene['Pos'], pdb_tool_gene["pACC"].fillna(0), 
                                label="pACC", zorder=3, color=sns.color_palette("tab10")[4], lw=0.5)      
                axes[ax].set_ylabel('Solvent\naccessibility', fontsize=13.5, rotation=0, va='center')
                axes[ax].yaxis.set_label_coords(-0.06, 0.5)
                axes[ax].set_ylim(-10, 110)
            
            # Plot stability change
            # ---------------------
            if "ddg" in annotations:
                ax = annotations.index("ddg")
                
                max_value, min_value = pos_result_gene["DDG"].max(), pos_result_gene["DDG"].min()
                # axes[ax+6].fill_between(pos_result_gene['Pos'], min_value, max_value, where=((pos_result_gene["C"] == 0) | (pos_result_gene["C"] == 2)), 
                #                         color='#ffd8b1', alpha=0.6)
                if sum(pos_result_gene['C'] == 1) > 0:
                    axes[ax].fill_between(pos_result_gene['Pos'], min_value, max_value, where=(pos_result_gene['C'] == 1), 
                                            color='white', lw=2)
                    axes[ax].fill_between(pos_result_gene['Pos'], min_value, max_value, where=(pos_result_gene['C'] == 1), 
                                            color='skyblue', alpha=0.4, lw=2)
                axes[ax].fill_between(pos_result_gene['Pos'], 0, pos_result_gene["DDG"], zorder=1,             
                                        color="white")     
                axes[ax].fill_between(pos_result_gene['Pos'], 0, pos_result_gene["DDG"], zorder=1,             
                                        color=sns.color_palette("pastel")[4], alpha=0.6)      
                axes[ax].plot(pos_result_gene['Pos'], pos_result_gene["DDG"], 
                                label="Stability change", zorder=2, color=sns.color_palette("tab10")[4], lw=0.5)    
                axes[ax].set_ylabel('ΔΔG (kcal/mol)', fontsize=13.5, rotation=0, va='center')
                axes[ax].yaxis.set_label_coords(-0.06, 0.5)
            
            # PTM
            # --------------
            if "ptm" in annotations:
                ax = annotations.index("ptm")

                ptm_gene = uni_feat_gene[uni_feat_gene["Type"] == "PTM"]
                ptm_names = ptm_gene["Description"].unique()
                sb_width = 0.5
                max_value = (len(ptm_names) * sb_width) - 0.2
                min_value = - 0.3
            
                # axes[ax].fill_between(pos_result_gene['Pos'], min_value, max_value, where=((pos_result_gene["C"] == 0) | (pos_result_gene["C"] == 2)), 
                #                 color='#ffd8b1', alpha=0.6, label='Mutated not *')
                axes[ax].fill_between(pos_result_gene['Pos'], min_value, max_value, where=(pos_result_gene['C'] == 1), 
                                        color='white', lw=2)
                axes[ax].fill_between(pos_result_gene['Pos'], min_value, max_value, where=(pos_result_gene['C'] == 1), 
                                        color='skyblue', alpha=0.4, label='Mutated *', lw=2)
            
                for n, name in enumerate(ptm_names):
                    c = sns.color_palette("tab10")[n]
                    ptm = ptm_gene[ptm_gene["Description"] == name]
                    ptm_pos = ptm.Begin.values
                    axes[ax].scatter(ptm_pos, np.repeat(n*sb_width, len(ptm_pos)), label=name, alpha=0.7, color=c) #label=name
                    axes[ax].hlines(y=n*sb_width, xmin=0, xmax=gene_len, linewidth=1, color='lightgray', alpha=0.7, zorder=0)
            
                axes[ax].set_ylim(min_value, max_value)
                y_ticks_positions = sb_width * np.arange(len(ptm_names))
                axes[ax].set_yticks(y_ticks_positions)
                axes[ax].set_yticklabels(ptm_names)
                axes[ax].set_ylabel(' PTM            ', fontsize=13.5, rotation=0, va='center')
            
            # SITES
            # --------------
            if "site" in annotations:
                ax = annotations.index("site")
            
                site_gene = uni_feat_gene[uni_feat_gene["Type"] == "SITE"]
                site_names = site_gene["Description"].unique()
                sb_width = 0.5
                max_value = (len(site_names) * sb_width) - 0.2
                min_value = - 0.3
                
                # axes[ax+8].fill_between(pos_result_gene['Pos'], min_value, max_value, where=((pos_result_gene["C"] == 0) | (pos_result_gene["C"] == 2)), 
                #                 color='#ffd8b1', alpha=0.6, label='Mutated not *')
                axes[ax].fill_between(pos_result_gene['Pos'], min_value, max_value, where=(pos_result_gene['C'] == 1), 
                                        color='white', lw=2)
                axes[ax].fill_between(pos_result_gene['Pos'], min_value, max_value, where=(pos_result_gene['C'] == 1), 
                                        color='skyblue', alpha=0.4, label='Mutated *', lw=2)
                
                for n, name in enumerate(site_names):
                    c = sns.color_palette("tab10")[n]
                    site_df = site_gene[site_gene["Description"] == name]
                    site_pos = get_site_pos(site_df)
                    axes[ax].scatter(site_pos, np.repeat(n*sb_width, len(site_pos)), label=name, alpha=0.7, color=c)
                    axes[ax].hlines(y=n*sb_width, xmin=0, xmax=gene_len, linewidth=1, color='lightgray', alpha=0.7, zorder=0)
                
                axes[ax].set_ylim(min_value, max_value)
                y_ticks_positions = sb_width * np.arange(len(site_names))
                axes[ax].set_yticks(y_ticks_positions)
                axes[ax].set_yticklabels(site_names)
                axes[ax].set_ylabel('Site           ', fontsize=13.5, rotation=0, va='center')
            
            # Clusters label
            # --------------
            if "clusters" in annotations:
                ax = annotations.index("clusters")

                clusters_label = pos_result_gene.Clump.dropna().unique()
                palette = sns.color_palette(cc.glasbey, n_colors=len(clusters_label))
                for i, cluster in enumerate(clusters_label):
                    axes[ax].fill_between(pos_result_gene['Pos'], -0.5, 0.46, 
                                            where=((pos_result_gene['Clump'] == cluster) & (pos_result_gene['C'] == 1)),
                                            color=palette[i], lw=0.4) # alpha=0.6
                axes[ax].set_ylabel('Clumps', fontsize=13.5, rotation=0, va='center')
                axes[ax].set_yticks([])  
                axes[ax
                ].set_yticklabels([], fontsize=12)
                axes[ax].yaxis.set_label_coords(-0.06, 0.5)
            
            # Secondary structure
            # -------------------
            if "sse" in annotations:
                ax = annotations.index("sse")
            
                for i, sse in enumerate(['Helix', 'Ladder', 'Coil']):
                    c = 0+i
                    ya, yb = c-plot_pars["sse_fill_width"], c+plot_pars["sse_fill_width"]
                    axes[ax].fill_between(pdb_tool_gene["Pos"].values, ya, yb, where=(pdb_tool_gene["SSE"] == sse), 
                                    color=sns.color_palette("tab10")[7+i], label=sse)
                axes[ax].set_yticks([0, 1, 2])  
                axes[ax].set_yticklabels(['Helix', 'Ladder', 'Coil'], fontsize=10)
                axes[ax].set_ylabel('SSE', fontsize=13.5, rotation=0, va='center')
                axes[ax].yaxis.set_label_coords(-0.06, 0.5)
            
            # Pfam
            # ----
            if "pfam" in annotations:
                ax = annotations.index("pfam")
            
                pfam_gene = uni_feat_gene[(uni_feat_gene["Type"] == "DOMAIN") & (uni_feat_gene["Evidence"] == "Pfam")]
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
                    axes[ax].fill_between(range(start, end+1), -0.45, 0.45,  alpha=0.5, color=pfam_color_dict[name])
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
                        axes[ax].text(((start + end) / 2)+0.5, y, name, ha='center', va='center', fontsize=10, color="black")
                        added_pfam.append(name)
                axes[ax].set_yticks([])  
                axes[ax].set_yticklabels([], fontsize=12)
                axes[ax].set_ylabel('Pfam', fontsize=13.5, rotation=0, va='center')
                axes[ax].set_ylim(-0.5, 0.5)  
                axes[ax].yaxis.set_label_coords(-0.06, 0.5)
            
            # Prosite
            # -------
            if "prosite" in annotations:
                ax = annotations.index("prosite")
            
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
                    axes[ax].fill_between(range(start, end+1), -0.45, 0.45,  alpha=0.5, color=prosite_color_dict[name])
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
                        axes[ax].text(((start + end) / 2)+0.5, y, name, ha='center', va='center', fontsize=10, color="black")
                        added_prosite.append(name)
                axes[ax].set_yticks([])  
                axes[ax].set_yticklabels([], fontsize=12)
                axes[ax].set_ylabel('Prosite', fontsize=13.5, rotation=0, va='center')
                axes[ax].set_ylim(-0.5, 0.5)  
                axes[ax].yaxis.set_label_coords(-0.06, 0.5)
            
            # Membrane
            # --------
            if "membrane" in annotations:
                ax = annotations.index("membrane")
            
                membrane_gene = uni_feat_gene[(uni_feat_gene["Type"] == "MEMBRANE")]
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
                    axes[ax].fill_between(range(start, end+1), -0.45, 0.45,  alpha=0.5, color=membrane_color_dict[name])
                    if name not in added_membrane:
                        y = -0.04
                        axes[ax].text(((start + end) / 2)+0.5, y, name, ha='center', va='center', fontsize=10, color="black")
                        added_membrane.append(name)
                axes[ax].set_yticks([])  
                axes[ax].set_yticklabels([], fontsize=12)
                axes[ax].set_ylabel('Membrane', fontsize=13.5, rotation=0, va='center')
                axes[ax].set_ylim(-0.5, 0.5)  
                axes[ax].yaxis.set_label_coords(-0.06, 0.5)
            
            # Motifs
            # ------
            if "motif" in annotations:
                ax = annotations.index("motif")
            
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
                    axes[ax].fill_between(range(start, end+1), -0.45, 0.45,  alpha=0.5, color=motif_color_dict[name])
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
                        axes[ax].text(((start + end) / 2)+0.5, y, name, ha='center', va='center', fontsize=10, color="black")
                        added_motif.append(name)
                axes[ax].set_yticks([])  
                axes[ax].set_yticklabels([], fontsize=12)
                axes[ax].set_ylabel('Motif', fontsize=13.5, rotation=0, va='center')
                axes[ax].set_ylim(-0.5, 0.5) 
                axes[ax].yaxis.set_label_coords(-0.06, 0.5)
                
            axes[len(axes)-1].set_xlabel(None)
            
            # Save
            # ====
            if title:
                fig.suptitle(f'{title}\n{gene} - {uni_id}', fontsize=16)
            else:
                fig.suptitle(f'{gene} - {uni_id}', fontsize=16)
            filename = f"{cohort}.genes_plot_{j+1}.{gene}_{uni_id}.png"
            output_path = os.path.join(output_dir, filename)
            htop = 0.947
            if title:
                htop -= 0.018
            plt.subplots_adjust(top=htop) 

            if save_plot:
                plt.savefig(output_path, dpi=300, bbox_inches='tight')
                logger.debug(f"Saved {output_path}")
            if show_plot: 
                plt.show()
            plt.close()

            # Store annotated result
            pos_result_gene = get_enriched_result(pos_result_gene, 
                                                  disorder_gene, 
                                                  pdb_tool_gene, 
                                                  seq_df)
            annotated_result_lst.append(pos_result_gene)
            uni_feat_result_lst.append(uni_feat_gene)

    pos_result_annotated = pd.concat(annotated_result_lst)
    feat_processed = pd.concat(uni_feat_result_lst)   
        
    return pos_result_annotated, feat_processed


# Comparative plots
# =================

def comparative_plots(shared_genes,
                    pos_result_1, 
                    maf_1,
                    maf_nonmiss_1,
                    miss_prob_dict_1,
                    cohort_1,
                    pos_result_2,
                    maf_2,
                    maf_nonmiss_2,
                    miss_prob_dict_2,
                    cohort_2,
                    seq_df,
                    output_dir,
                    annotations_dir,
                    disorder,
                    uniprot_feat,
                    pdb_tool,
                    plot_pars,
                    save_plot=True,
                    show_plot=False):   
    """
    Generate plot to compare each gene that are processed in both 3D-clustering analysis.
    """

    warnings.filterwarnings("ignore", category=UserWarning)
    
    for j, gene in enumerate(shared_genes):
        
        logger.debug(f"Generating comparative plots for {len(shared_genes)} genes..")
      
        # Load and parse
        # ==============
        
        uni_id = seq_df[seq_df["Gene"] == gene].Uniprot_ID.values[0]
        af_f = seq_df[seq_df["Gene"] == gene].F.values[0]
        gene_len = len(seq_df[seq_df["Gene"] == gene].Seq.values[0])
        maf_gene_1 = maf_1[maf_1["Gene"] == gene]
        maf_gene_2 = maf_2[maf_2["Gene"] == gene]
        
        # Parse
        pos_result_gene_1 = pos_result_1[pos_result_1["Gene"] == gene].sort_values("Pos").reset_index(drop=True)
        pos_result_gene_2 = pos_result_2[pos_result_2["Gene"] == gene].sort_values("Pos").reset_index(drop=True)

        if len(pos_result_gene_1) > 0 and len(pos_result_gene_2) > 0:
            pos_result_gene_1 = pos_result_gene_1[["Pos", "Mut_in_res", "Mut_in_vol", 
                                                  "Score_obs_sim", "C", "C_ext", 
                                                  "pval", "Clump", "PAE_vol"]]
            pos_result_gene_2 = pos_result_gene_2[["Pos", "Mut_in_res", "Mut_in_vol", 
                                                   "Score_obs_sim", "C", "C_ext", 
                                                   "pval", "Clump", "PAE_vol"]]
            pos_result_gene_1, max_mut_1 = parse_pos_result_for_genes_plot(pos_result_gene_1)
            pos_result_gene_2, max_mut_2 = parse_pos_result_for_genes_plot(pos_result_gene_2)
            
            # Counts
            mut_count_1, mut_count_nonmiss_1 = get_count_for_genes_plot(maf_gene_1, 
                                                                        maf_nonmiss_1, 
                                                                        gene, 
                                                                        non_missense_count="nonmiss_count" in plot_pars["h_ratios"])
            mut_count_2, mut_count_nonmiss_2 = get_count_for_genes_plot(maf_gene_2, 
                                                                        maf_nonmiss_2, 
                                                                        gene, 
                                                                        non_missense_count="nonmiss_count" in plot_pars["h_ratios"])
            
            # Get prob vec
            prob_vec_1 = np.array(miss_prob_dict_1[f"{uni_id}-F{af_f}"])  
            prob_vec_2 = np.array(miss_prob_dict_2[f"{uni_id}-F{af_f}"])    
        
            # Get per-pos score and normalize score
            pos_result_gene_1, score_vec_1, score_norm_vec_1 = get_score_for_genes_plot(pos_result_gene_1, 
                                                                                        mut_count_1, 
                                                                                        prob_vec_1)
            pos_result_gene_2, score_vec_2, score_norm_vec_2 = get_score_for_genes_plot(pos_result_gene_2, 
                                                                                        mut_count_2, 
                                                                                        prob_vec_2)

            # Get annotations
            pos_result_gene_1, disorder_gene, pdb_tool_gene, uni_feat_gene = get_id_annotations(uni_id, 
                                                                                                pos_result_gene_1, 
                                                                                                maf_gene_1, 
                                                                                                annotations_dir, 
                                                                                                disorder, 
                                                                                                pdb_tool, 
                                                                                                uniprot_feat)
            pos_result_gene_2, _, _, _ = get_id_annotations(uni_id, 
                                                            pos_result_gene_2, 
                                                            maf_gene_2, 
                                                            annotations_dir, 
                                                            disorder, 
                                                            pdb_tool, 
                                                            uniprot_feat)

            # Pos result for background filling
            pos_result_gene_shared = pos_result_gene_1.copy()
            pos_result_gene_shared["C"] = np.nan
            pos_result_gene_shared["C_A"] = pos_result_gene_1["C"]
            pos_result_gene_shared["C_B"] = pos_result_gene_2["C"]
            pos_result_gene_shared["C"] = pos_result_gene_shared.apply(lambda x: 
                                                                       "A" if x.C_A == 1 and x.C_B != 1 else 
                                                                       "B" if x.C_A != 1 and x.C_B == 1 else 
                                                                       "AB" if x.C_A == 1 and x.C_B == 1 else np.nan, axis=1)
            
            # Generate plot
            # ============= 
                
            if not maf_nonmiss_1 or not maf_nonmiss_2:
                maf_nonmiss = None
            else:
                maf_nonmiss = maf_nonmiss_1
            h_ratios, near_pfam, near_prosite, near_motif = get_gene_arg(pd.concat((pos_result_gene_1, pos_result_gene_2)), 
                                                                         plot_pars, 
                                                                         uni_feat_gene, 
                                                                         maf_nonmiss)
            annotations = list(h_ratios.keys())
            
            fig, axes = plt.subplots(len(h_ratios), 1, 
                                     figsize=plot_pars["figsize"], 
                                     sharex=True, 
                                     gridspec_kw={'hspace': 0.1, 
                                                  'height_ratios': h_ratios.values()})
                
                
            # Plot for Non-missense mut track            ## TO DO: Enable not mirror for non-missense
            # -------------------------------
            if "nonmiss_count" in annotations:
                ax = annotations.index("nonmiss_count")
            
                if len(mut_count_nonmiss_1.Consequence.unique()) > 6:
                    ncol = 3
                else:
                    ncol = 2
                i = 0
                axes[ax].vlines(mut_count_nonmiss_1["Pos"], ymin=0, ymax=mut_count_nonmiss_1["Count"], 
                                color="gray", lw=0.7, zorder=0, alpha=0.5) # To cover the overlapping needle top part
                axes[ax].scatter(mut_count_nonmiss_1["Pos"], mut_count_nonmiss_1["Count"], color='white', zorder=4, lw=plot_pars["s_lw"]) 
                for cnsq in mut_count_nonmiss_1.Consequence.unique():
                    count_cnsq = mut_count_nonmiss_1[mut_count_nonmiss_1["Consequence"] == cnsq]
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
                                    color=color, zorder=order, alpha=0.7, lw=plot_pars["s_lw"], ec="black")              # ec="black",
                axes[ax].legend(fontsize=11.5, ncol=ncol, framealpha=0.75)
                axes[ax].set_ylabel('Non\nmissense\nmutations', fontsize=13.5, rotation=0, va='center')
                ymargin = max(max(mut_count_nonmiss_1["Count"]), max(mut_count_nonmiss_2["Count"])) * 0.1
                axes[ax].set_ylim(-(max(mut_count_nonmiss_1["Count"])-ymargin), max(mut_count_nonmiss_1["Count"])+ymargin)


            # Plot for Missense mut track
            # ---------------------------
            mut_pos_1 = pos_result_gene_1[pos_result_gene_1["Mut_in_res"] > 0].Pos.values
            mut_res_pos_1 = pos_result_gene_1[pos_result_gene_1["Mut_in_res"] > 0].Mut_in_res.values
            mut_pos_2 = pos_result_gene_2[pos_result_gene_2["Mut_in_res"] > 0].Pos.values
            mut_res_pos_2 = pos_result_gene_2[pos_result_gene_2["Mut_in_res"] > 0].Mut_in_res.values

            if plot_pars["count_mirror"]:
                if "miss_count" in annotations:
                    ax = annotations.index("miss_count")
        
                    axes[ax].hlines(0, xmin=0, xmax=gene_len, color="gray", lw=0.6, zorder=1)
                    axes[ax].vlines(mut_count_1["Pos"], ymin=0, ymax=mut_count_1["Count"], color="gray", lw=0.7, zorder=1, alpha=0.5)            # A
                    axes[ax].vlines(mut_count_2["Pos"], ymin=-mut_count_2["Count"], ymax=0, color="gray", lw=0.7, zorder=1, alpha=0.5)       # B

                    axes[ax].fill_between(pos_result_gene_1['Pos'], 0, 0, where=(pos_result_gene_1['C'] == "NA"), 
                                            color=sns.color_palette("pastel")[2], alpha=0.4, label='Position in cluster A', zorder=0, lw=2) # Just for the legend
                    axes[ax].fill_between(pos_result_gene_1['Pos'], 0, 0, where=(pos_result_gene_1['C'] == "NA"), 
                                            color=sns.color_palette("pastel")[3], alpha=0.4, label='Position in cluster B', zorder=0, lw=2) # Just for the legend
                    axes[ax].fill_between(pos_result_gene_shared['Pos'], -max_mut_2, max_mut_1, 
                                            where=(pos_result_gene_shared['C'] == "A") | (pos_result_gene_shared['C'] == "B") | (pos_result_gene_shared['C'] == "AB"), 
                                            color='skyblue', alpha=0.4, label='Position in cluster A or B', zorder=0, lw=2)
                    
                    axes[ax].scatter(mut_pos_1, mut_res_pos_1, color='white', zorder=3, lw=plot_pars["s_lw"], ec="white")               # A
                    axes[ax].scatter(mut_pos_2, -mut_res_pos_2, color='white', zorder=3, lw=plot_pars["s_lw"], ec="white")          # B
                    
                    axes[ax].scatter(mut_pos_1, mut_res_pos_1, color="C2", zorder=4, alpha=0.6,                   # A
                                    lw=plot_pars["s_lw"], ec="black", s=60, label='Cohort A')    
                    axes[ax].scatter(mut_pos_2, -mut_res_pos_2, color="tomato", zorder=4, alpha=0.6,          # B
                                    lw=plot_pars["s_lw"], ec="black", s=60, label='Cohort B')    
        

                    legend = axes[ax].legend(fontsize=11.5, ncol=2, framealpha=0.75, bbox_to_anchor=(0.95, 2.3), 
                                             borderaxespad=0., loc='upper right')
                    legend.set_title("Global legend")
                    legend.get_title().set_fontsize(12)
                    
                    axes[ax].set_ylabel('Missense\nmutations', fontsize=13.5, rotation=0, va='center') 
                    axes[ax].yaxis.set_label_coords(-0.06, 0.5)
                    ymargin = max(max(mut_res_pos_2), max(mut_res_pos_1)) * 0.1
                    axes[ax].set_ylim(-max(mut_res_pos_2)-ymargin, max(mut_res_pos_1)+ymargin)

                    tick_labels = [f'{abs(label):.4g}' for label in axes[ax].get_yticks()]
                    axes[ax].set_yticklabels(tick_labels)

            else:
                if "miss_count" in annotations and "miss_count_2" in annotations:
                       
                    # A
                    ax = annotations.index("miss_count")
                    axes[ax].vlines(mut_count_1["Pos"], ymin=0, ymax=mut_count_1["Count"], color="gray", lw=0.7, zorder=1, alpha=0.5)           
                    axes[ax].fill_between(pos_result_gene_1['Pos'], 0, max(mut_res_pos_1), where=(pos_result_gene_1['C'] == 1), 
                                            color=sns.color_palette("pastel")[2], alpha=0.4, label='Position in cluster A', zorder=0, lw=2)
                    axes[ax].fill_between(pos_result_gene_1['Pos'], 0, 0, where=(pos_result_gene_1['C'] == "NA"), 
                                            color=sns.color_palette("pastel")[3], alpha=0.4, label='Position in cluster B', zorder=0, lw=2) # Just for the legend
                    axes[ax].fill_between(pos_result_gene_1['Pos'], 0, 0, where=(pos_result_gene_1['C'] == "NA"), 
                                            color="skyblue", alpha=0.4, label='Position in cluster A or B', zorder=0, lw=2) # Just for the legend
                    
                    axes[ax].scatter(mut_pos_1, mut_res_pos_1, color='white', zorder=3, lw=plot_pars["s_lw"], ec="white")            
                    axes[ax].scatter(mut_pos_1, mut_res_pos_1, color="C2", zorder=4, alpha=0.6,                  
                                    lw=plot_pars["s_lw"], ec="black", s=60, label='Cohort A')
                    axes[ax].scatter(-20, -20, color="tomato", zorder=4, alpha=0.6,             # Just for the legend                           
                                    lw=plot_pars["s_lw"], ec="black", s=60, label='Cohort B')
                    
                    legend = axes[ax].legend(fontsize=11.5, ncol=2, framealpha=0.75, 
                                            bbox_to_anchor=(0.95, 2.45), borderaxespad=0., loc='upper right')
                    legend.set_title("Global legend")
                    legend.get_title().set_fontsize(12)

                    axes[ax].set_ylabel('Missense\nmutations A', fontsize=13.5, rotation=0, va='center') 
                    axes[ax].yaxis.set_label_coords(-0.06, 0.5)
                    ymargin = max(mut_res_pos_1) * 0.1
                    axes[ax].set_ylim(0-ymargin, max(mut_res_pos_1)+ymargin)

                    # B
                    ax = annotations.index("miss_count_2")
                    axes[ax].vlines(mut_count_2["Pos"], ymin=0, ymax=mut_count_2["Count"], color="gray", lw=0.7, zorder=1, alpha=0.5)      
                    axes[ax].fill_between(pos_result_gene_2['Pos'], 0, max(mut_res_pos_2), where=(pos_result_gene_2['C'] == 1), 
                                            color=sns.color_palette("pastel")[3], alpha=0.4, label='Position in cluster B', zorder=0, lw=2)
                    axes[ax].scatter(mut_pos_2, mut_res_pos_2, color='white', zorder=3, lw=plot_pars["s_lw"], ec="white")         
                    axes[ax].scatter(mut_pos_2, mut_res_pos_2, color="tomato", zorder=4, alpha=0.6,        
                                    lw=plot_pars["s_lw"], ec="black", s=60, label='Cohort B')  
                    axes[ax].set_ylabel('Missense\nmutations B', fontsize=13.5, rotation=0, va='center') 
                    axes[ax].yaxis.set_label_coords(-0.06, 0.5)
                    ymargin = max(mut_res_pos_2) * 0.1
                    axes[ax].set_ylim(0-ymargin, max(mut_res_pos_2)+ymargin)


            # Plot for Miss prob track
            # ------------------------
            if "miss_prob" in annotations:
                ax = annotations.index("miss_prob")

                if plot_pars["prob_mirror"]:
                    max_value = max(prob_vec_1)
                    min_value = -max(prob_vec_2)
                    prob_vec_2 = -np.array(prob_vec_2)
                else:
                    max_value = max(max(prob_vec_2), max(prob_vec_1))
                    min_value = 0

                axes[ax].fill_between(pos_result_gene_shared['Pos'], min_value, max_value, 
                                        where=(pos_result_gene_shared['C'] == "A") | (pos_result_gene_shared['C'] == "B") | (pos_result_gene_shared['C'] == "AB"), 
                                        color='skyblue', alpha=0.4, label='Position in cluster', zorder=0, lw=2)

                axes[ax].hlines(0, xmin=0, xmax=gene_len, color="gray", lw=0.6, zorder=1)
                axes[ax].plot(range(1, len(prob_vec_1)+1), prob_vec_1, label="Cohort A", zorder=3, color="C2", lw=1)                          
                axes[ax].plot(range(1, len(prob_vec_2)+1), prob_vec_2, label="Cohort B", zorder=3, 
                                color="tomato", lw=1)                          
                
                axes[ax].set_ylabel('Missense\nmut prob', fontsize=13.5, rotation=0, va='center')
                axes[ax].yaxis.set_label_coords(-0.06, 0.5)
                if plot_pars["prob_mirror"]:
                    tick_labels = [f'{abs(label):.4g}' for label in axes[ax].get_yticks()]
                    axes[ax].set_yticklabels(tick_labels)
                
            
            # Plot for Score track
            # --------------------
            if plot_pars["score_mirror"]:
                
                if "score" in annotations:
                    ax = annotations.index("score")
                    
                    max_value = np.max(score_vec_1)
                    min_value = -np.max(score_vec_2)
                    axes[ax].fill_between(pos_result_gene_shared['Pos'], min_value, max_value, 
                                            where=(pos_result_gene_shared['C'] == "A") | (pos_result_gene_shared['C'] == "B") | (pos_result_gene_shared['C'] == "AB"), 
                                            color='skyblue', alpha=0.4, label='Position in cluster', zorder=0, lw=2)
                    axes[ax].hlines(0, xmin=0, xmax=gene_len, color="gray", lw=0.7, zorder=1)
                    axes[ax].plot(range(1, len(prob_vec_1)+1), score_vec_1, label="Cohort A", zorder=2, color="C2", lw=1)                       
                    axes[ax].plot(range(1, len(prob_vec_2)+1), -np.array(score_vec_2), label="Cohort B", zorder=2, color="tomato", lw=1)  
                    
                    axes[ax].set_ylabel('Clustering\nscore\n(obs/sim)', fontsize=13.5, rotation=0, va='center')
                    axes[ax].yaxis.set_label_coords(-0.06, 0.5)
                    
                    tick_labels = [f'{abs(label):.4g}' for label in axes[ax].get_yticks()]
                    axes[ax].set_yticklabels(tick_labels)

            else:
                if "score" in annotations and "score_2" in annotations:
                    
                    # A
                    ax = annotations.index("score")
                    
                    max_value = np.max(score_vec_1)
                    axes[ax].fill_between(pos_result_gene_1['Pos'], 0, max_value, where=(pos_result_gene_1['C'] == 1), 
                                            color=sns.color_palette("pastel")[2], alpha=0.4, label='Position in cluster A', zorder=0, lw=2)
                    axes[ax].plot(range(1, len(prob_vec_1)+1), score_vec_1, label="Cohort A", zorder=2, color="C2", lw=1)                       
                    axes[ax].set_ylabel('Clustering\nscore A\n(obs/sim)', fontsize=13.5, rotation=0, va='center')
                    axes[ax].yaxis.set_label_coords(-0.06, 0.5)

                    # B
                    ax = annotations.index("score_2")
                    
                    max_value = np.max(score_vec_2)
                    axes[ax].fill_between(pos_result_gene_2['Pos'], 0, max_value, where=(pos_result_gene_2['C'] == 1), 
                                            color=sns.color_palette("pastel")[3], alpha=0.4, label='Position in cluster B', zorder=0, lw=2)
                    axes[ax].plot(range(1, len(prob_vec_2)+1), np.array(score_vec_2), label="Cohort B", zorder=2, color="tomato", lw=1)  
                    axes[ax].set_ylabel('Clustering\nscore B\n(obs/sim)', fontsize=13.5, rotation=0, va='center') 
                    axes[ax].yaxis.set_label_coords(-0.06, 0.5)


            # Clusters label A
            # ----------------
            
            # A
            if "clusters" in annotations: 
                ax = annotations.index("clusters") 

                clusters_label = pos_result_gene_1.Clump.dropna().unique()
                clusters_label_2 = pos_result_gene_2.Clump.dropna().unique()
                n_colors = max(len(clusters_label), len(clusters_label_2))
                palette = sns.color_palette(cc.glasbey, n_colors=n_colors)
                for i, cluster in enumerate(clusters_label):
                    axes[ax].fill_between(pos_result_gene_1['Pos'], -0.5, 0.46, 
                                            where=((pos_result_gene_1['Clump'] == cluster) & (pos_result_gene_1['C'] == 1)),
                                            color=palette[i], lw=0.4) # alpha=0.6
                axes[ax].set_ylabel('Clusters A                ', fontsize=13.5, rotation=0, va='center')
                axes[ax].set_yticks([])  
                axes[ax].yaxis.set_label_coords(-0.034, 0.5)

            # B
            if "clusters_2" in annotations: 
                ax = annotations.index("clusters_2") 
                
                clusters_label_2 = pos_result_gene_2.Clump.dropna().unique()
                for i, cluster in enumerate(clusters_label_2):
                    axes[ax].fill_between(pos_result_gene_2['Pos'], -0.5, 0.46, 
                                            where=((pos_result_gene_2['Clump'] == cluster) & (pos_result_gene_2['C'] == 1)),
                                            color=palette[i], lw=0.4) # alpha=0.6
                axes[ax].set_ylabel('Clusters B                ', fontsize=13.5, rotation=0, va='center')
                axes[ax].set_yticks([])  
                axes[ax].yaxis.set_label_coords(-0.034, 0.5)


            # Plot annotations
            # ================

            # Plot PAE
            # --------
            if "pae" in annotations: 
                ax = annotations.index("pae")
                
                max_value = np.max(pos_result_gene_1["PAE_vol"])
                axes[ax].fill_between(pos_result_gene_shared['Pos'], 0, max_value, 
                                where=(pos_result_gene_shared['C'] == "A") | (pos_result_gene_shared['C'] == "B") | (pos_result_gene_shared['C'] == "AB"), 
                                color='skyblue', alpha=0.4, label='Position in cluster', zorder=0, lw=2)
                axes[ax].fill_between(pos_result_gene_1["Pos"], 0, pos_result_gene_1["PAE_vol"].fillna(0), 
                                        zorder=2, color="white")    
                axes[ax].fill_between(pos_result_gene_1["Pos"], 0, pos_result_gene_1["PAE_vol"].fillna(0), 
                                        zorder=2, color=sns.color_palette("pastel")[4], alpha=0.6)    
                axes[ax].plot(pos_result_gene_1['Pos'], pos_result_gene_1["PAE_vol"].fillna(0),                                     
                                label="Confidence", zorder=3, color=sns.color_palette("tab10")[4], lw=0.5)
                axes[ax].set_ylabel('Predicted\naligned\nerror\n(Å)', fontsize=13.5, rotation=0, va='center')
                axes[ax].yaxis.set_label_coords(-0.06, 0.5)

                
            # Plot disorder
            # -------------
            if "disorder" in annotations: 
                ax = annotations.index("disorder")
                
                axes[ax].fill_between(pos_result_gene_shared['Pos'], 0, 100, 
                                        where=(pos_result_gene_shared['C'] == "A") | (pos_result_gene_shared['C'] == "B") | (pos_result_gene_shared['C'] == "AB"), 
                                        color='skyblue', alpha=0.4, label='Position in cluster', zorder=0, lw=2)

                # ## Comment out to use AF color palette
                
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
                #     axes[ax].fill_between(disorder_x, 0, disorder_y, where=(condition),       
                #                             zorder=2, color="white")   
                #     axes[ax].fill_between(disorder_x, 0, disorder_y, where=(condition),   
                #                             zorder=3, facecolor=color, alpha=0.8)  

                axes[ax].fill_between(disorder_gene["Pos"], 0, disorder_gene["Confidence"].fillna(0),                  
                                        zorder=2, color="white")
                axes[ax].fill_between(disorder_gene["Pos"], 0, disorder_gene["Confidence"].fillna(0),                  
                                        zorder=2, color=sns.color_palette("pastel")[4], alpha=0.6)
            
                
                axes[ax].plot(disorder_gene["Pos"], disorder_gene["Confidence"], 
                                label="Confidence", zorder=3, color=sns.color_palette("tab10")[4], lw=0.5)    
                axes[ax].set_ylabel('pLDDT\n(disorder)', fontsize=13.5, rotation=0, va='center')
                axes[ax].yaxis.set_label_coords(-0.06, 0.5)
                axes[ax].set_ylim(-10, 110)
                    

            # Plot pACC
            # ---------
            if "pacc" in annotations: 
                ax = annotations.index("pacc")

                axes[ax].fill_between(pos_result_gene_shared['Pos'], 0, 100, 
                                        where=(pos_result_gene_shared['C'] == "A") | (pos_result_gene_shared['C'] == "B") | (pos_result_gene_shared['C'] == "AB"), 
                                        color='skyblue', alpha=0.4, label='Position in cluster', zorder=0, lw=2) 
                axes[ax].fill_between(pdb_tool_gene["Pos"], 0, pdb_tool_gene["pACC"].fillna(0),                  
                                        zorder=2, color="white")
                axes[ax].fill_between(pdb_tool_gene["Pos"], 0, pdb_tool_gene["pACC"].fillna(0),                  
                                        zorder=2, color=sns.color_palette("pastel")[4], alpha=0.6)
                axes[ax].plot(pdb_tool_gene['Pos'], pdb_tool_gene["pACC"].fillna(0), 
                                label="pACC", zorder=3, color=sns.color_palette("tab10")[4], lw=0.5)      
                axes[ax].set_ylabel('Solvent\naccessibility', fontsize=13.5, rotation=0, va='center')
                axes[ax].yaxis.set_label_coords(-0.06, 0.5)
                axes[ax].set_ylim(-10, 110)


            # Plot stability change A
            # -----------------------
            if "ddg" in annotations: 
                ax = annotations.index("ddg")

                max_value, min_value = pos_result_gene_1["DDG"].max(), pos_result_gene_1["DDG"].min()

                axes[ax].fill_between(pos_result_gene_1['Pos'], min_value, max_value, where=(pos_result_gene_1['C'] == 1), 
                                    color=sns.color_palette("pastel")[2], alpha=0.4, label='Position in cluster A', zorder=0, lw=2)

                axes[ax].fill_between(pos_result_gene_1['Pos'], 0, pos_result_gene_1["DDG"], zorder=1,             
                                        color="white")     
                axes[ax].fill_between(pos_result_gene_1['Pos'], 0, pos_result_gene_1["DDG"], zorder=1,             
                                        color=sns.color_palette("pastel")[4], alpha=0.6)      
                axes[ax].plot(pos_result_gene_1['Pos'], pos_result_gene_1["DDG"], 
                                label="Stability change", zorder=2, color=sns.color_palette("tab10")[4], lw=0.5)    
                axes[ax].set_ylabel('ΔΔG A\n(kcal/mol)', fontsize=13.5, rotation=0, va='center')
                axes[ax].yaxis.set_label_coords(-0.06, 0.5)


            # Plot stability change B
            # -----------------------
            if "ddg_2" in annotations: 
                ax = annotations.index("ddg_2")

                max_value, min_value = pos_result_gene_2["DDG"].max(), pos_result_gene_1["DDG"].min()

                axes[ax].fill_between(pos_result_gene_2['Pos'], min_value, max_value, where=(pos_result_gene_2['C'] == 1), 
                                        color=sns.color_palette("pastel")[3], alpha=0.4, label='Position in cluster B', zorder=0, lw=2)

                axes[ax].fill_between(pos_result_gene_2['Pos'], 0, pos_result_gene_2["DDG"], zorder=1,             
                                        color="white")     
                axes[ax].fill_between(pos_result_gene_2['Pos'], 0, pos_result_gene_2["DDG"], zorder=1,             
                                        color=sns.color_palette("pastel")[4], alpha=0.6)      
                axes[ax].plot(pos_result_gene_2['Pos'], pos_result_gene_2["DDG"], 
                                label="Stability change", zorder=2, color=sns.color_palette("tab10")[4], lw=0.5)    
                axes[ax].set_ylabel('ΔΔG B\n(kcal/mol)', fontsize=13.5, rotation=0, va='center')
                axes[ax].yaxis.set_label_coords(-0.06, 0.5)


            # PTM
            # ---
            if "ptm" in annotations: 
                ax = annotations.index("ptm")
        
                ptm_gene = uni_feat_gene[uni_feat_gene["Type"] == "PTM"]
                ptm_names = ptm_gene["Description"].unique()
                sb_width = 0.5
                max_value = (len(ptm_names) * sb_width) - 0.2
                min_value = - 0.3

                axes[ax].fill_between(pos_result_gene_shared['Pos'], min_value, max_value, 
                                        where=(pos_result_gene_shared['C'] == "A") | (pos_result_gene_shared['C'] == "B") | (pos_result_gene_shared['C'] == "AB"), 
                                        color='skyblue', alpha=0.4, label='Position in cluster', zorder=0, lw=2)

                for n, name in enumerate(ptm_names):
                    c = sns.color_palette("tab10")[n]
                    ptm = ptm_gene[ptm_gene["Description"] == name]
                    ptm_pos = ptm.Begin.values
                    axes[ax].scatter(ptm_pos, np.repeat(n*sb_width, len(ptm_pos)), label=name, alpha=0.7, color=c) #label=name
                    axes[ax].hlines(y=n*sb_width, xmin=0, xmax=gene_len, linewidth=1, color='lightgray', alpha=0.7, zorder=0)
            
                axes[ax].set_ylim(min_value, max_value)
                y_ticks_positions = sb_width * np.arange(len(ptm_names))
                axes[ax].set_yticks(y_ticks_positions)
                axes[ax].set_yticklabels(ptm_names)
                axes[ax].set_ylabel(' PTM            ', fontsize=13.5, rotation=0, va='center')


            # SITES
            # --------------
            if "site" in annotations: 
                ax = annotations.index("site") 

                site_gene = uni_feat_gene[uni_feat_gene["Type"] == "SITE"]
                site_names = site_gene["Description"].unique()
                sb_width = 0.5
                max_value = (len(site_names) * sb_width) - 0.2
                min_value = - 0.3

                axes[ax].fill_between(pos_result_gene_shared['Pos'], min_value, max_value, 
                                        where=(pos_result_gene_shared['C'] == "A") | (pos_result_gene_shared['C'] == "B") | (pos_result_gene_shared['C'] == "AB"), 
                                        color='skyblue', alpha=0.4, label='Position in cluster', zorder=0, lw=2)
                
                for n, name in enumerate(site_names):
                    c = sns.color_palette("tab10")[n]
                    site = site_gene[site_gene["Description"] == name]
                    site_pos = site.Begin.values
                    axes[ax].scatter(site_pos, np.repeat(n*sb_width, len(site_pos)), label=name, alpha=0.7, color=c) #label=name
                    axes[ax].hlines(y=n*sb_width, xmin=0, xmax=gene_len, linewidth=1, color='lightgray', alpha=0.7, zorder=0)
                
                axes[ax].set_ylim(min_value, max_value)
                y_ticks_positions = sb_width * np.arange(len(site_names))
                axes[ax].set_yticks(y_ticks_positions)
                axes[ax].set_yticklabels(site_names)
                axes[ax].set_ylabel('Site           ', fontsize=13.5, rotation=0, va='center')

                
            # Secondary structure
            # -------------------
            if "sse" in annotations: 
                ax = annotations.index("sse") 

                for i, sse in enumerate(['Helix', 'Ladder', 'Coil']):
                    c = 0+i
                    ya, yb = c-plot_pars["sse_fill_width"], c+plot_pars["sse_fill_width"]
                    axes[ax].fill_between(pdb_tool_gene["Pos"].values, ya, yb, where=(pdb_tool_gene["SSE"] == sse), 
                                    color=sns.color_palette("tab10")[7+i], label=sse)
                axes[ax].set_yticks([0, 1, 2])  
                axes[ax].set_yticklabels(['Helix', 'Ladder', 'Coil'], fontsize=10)
                axes[ax].set_ylabel('SSE       ', fontsize=13.5, rotation=0, va='center')
                axes[ax].yaxis.set_label_coords(-0.051, 0.5)


            # Pfam
            # ----
            if "pfam" in annotations: 
                ax = annotations.index("pfam") 
                
                pfam_gene = uni_feat_gene[(uni_feat_gene["Type"] == "DOMAIN") & (uni_feat_gene["Evidence"] == "Pfam")]
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
                    axes[ax].fill_between(range(start, end+1), -0.45, 0.45,  alpha=0.5, color=pfam_color_dict[name])
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
                        axes[ax].text(((start + end) / 2)+0.5, y, name, ha='center', va='center', fontsize=10, color="black")
                        added_pfam.append(name)
                axes[ax].set_yticks([])  
                axes[ax].set_yticklabels([], fontsize=12)
                axes[ax].set_ylabel('Pfam        ', fontsize=13.5, rotation=0, va='center')
                axes[ax].set_ylim(-0.5, 0.5)  
                axes[ax].yaxis.set_label_coords(-0.051, 0.5)


            # Prosite
            # -------
            if "prosite" in annotations: 
                ax = annotations.index("prosite") 
           
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
                    axes[ax].fill_between(range(start, end+1), -0.45, 0.45,  alpha=0.5, color=prosite_color_dict[name])
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
                        axes[ax].text(((start + end) / 2)+0.5, y, name, ha='center', va='center', fontsize=10, color="black")
                        added_prosite.append(name)
                axes[ax].set_yticks([])  
                axes[ax].set_yticklabels([], fontsize=12)
                axes[ax].set_ylabel('Prosite           ', fontsize=13.5, rotation=0, va='center')
                axes[ax].set_ylim(-0.5, 0.5)  
                axes[ax].yaxis.set_label_coords(-0.051, 0.5)


            # Membrane
            # --------
            if "membrane" in annotations: 
                ax = annotations.index("membrane") 

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
                    axes[ax].fill_between(range(start, end+1), -0.45, 0.45,  alpha=0.5, color=membrane_color_dict[name])
                    if name not in added_membrane:
                        y = -0.04
                        axes[ax].text(((start + end) / 2)+0.5, y, name, ha='center', va='center', fontsize=10, color="black")
                        added_membrane.append(name)
                axes[ax].set_yticks([])  
                axes[ax].set_yticklabels([], fontsize=12)
                axes[ax].set_ylabel('Membrane        ', fontsize=13.5, rotation=0, va='center')
                axes[ax].set_ylim(-0.5, 0.5)  
                axes[ax].yaxis.set_label_coords(-0.051, 0.5)
                    

            # Motifs
            # ------
            if "motif" in annotations: 
                ax = annotations.index("motif") 
                
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
                    axes[ax].fill_between(range(start, end+1), -0.45, 0.45,  alpha=0.5, color=motif_color_dict[name])
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
                        axes[ax].text(((start + end) / 2)+0.5, y, name, ha='center', va='center', fontsize=10, color="black")
                        added_motif.append(name)
                axes[ax].set_yticks([])  
                axes[ax].set_yticklabels([], fontsize=12)
                axes[ax].set_ylabel('Motif        ', fontsize=13.5, rotation=0, va='center')
                axes[ax].set_ylim(-0.5, 0.5) 
                axes[ax].yaxis.set_label_coords(-0.051, 0.5)


            # Save
            # ====
            
            fig.suptitle(f'{cohort_1} (A) - {cohort_2} (B)\n\n{gene} ({uni_id})', fontsize=16)
            if save_plot:
                filename = f"{cohort_1}.{cohort_2}.comp_plot_{j+1}.{gene}_{uni_id}.png"
                output_path = os.path.join(output_dir, filename)
                plt.subplots_adjust(top=0.9) 
                plt.savefig(output_path, dpi=300, bbox_inches='tight')
                logger.debug(f"Saved {output_path}")
            if show_plot: 
                plt.show()
            plt.close()
        
        else:
            logger.warning("Nothing to plot!")  
            
            
            
# Associations plots
# ==================

def expand_uniprot_feat_rows(df):
    """
    Convert Uniprot features df from ranges of positions to individual positions.
    """
    
    expanded_list = []
    for _, row in df.iterrows():
        begin, end = int(row['Begin']), int(row['End'])
        expanded_data = [{'Uniprot_ID': row['Uniprot_ID'], 'Type': row['Type'], 'Pos': pos, 'Description': row['Description']} 
                         for pos in range(begin, end + 1)]
        expanded_list.extend(expanded_data)
    expanded_df = pd.DataFrame(expanded_list)
    
    return expanded_df

    
def get_dummy(df, pos, annot):

    return int(annot in df[df["Pos"] == pos].Type.values)


def get_dummies_annot(df, col):

    pos = pd.Series(df.Pos.unique(), name="Pos")
    lst_results = [pos]
    for annot in df[col].unique():
        series = pd.Series(df.Pos.unique()).apply(lambda x: get_dummy(df, pos=x, annot=annot))
        series.name = annot.capitalize()
        lst_results.append(series)

    return pd.concat(lst_results, axis=1)


def get_uni_feat_for_odds(uni_id, uni_feat_df):
    """
    Preprocess Uniprot features df for logistic regression.
    """

    uni_feat_df = uni_feat_df.copy()
    uni_feat_df = uni_feat_df[uni_feat_df["Uniprot_ID"] == uni_id]
    uni_feat_df = get_dummies_annot(uni_feat_df, col="Type")
    uni_feat_df.insert(0, "Uniprot_ID", uni_id)

    return uni_feat_df
    

def uni_log_reg(df, labels):
    """
    Univariate logistic regression analysis.
    """

    df = df.copy()
    results = {}
    
    df = df.drop(columns=[col for col in df.columns if df[col].nunique() == 1])
    columns = df.columns
    
    # Keep tracks of NA in ddG before standardizing
    if "ΔΔG" in columns:
        ddg_ix = ~df["ΔΔG"].isna().values
    df = df.fillna(0)

    scaler = StandardScaler()
    df = scaler.fit_transform(df) 
    
    for i, col in enumerate(columns):
        
        # Drop NA only in since it is the only annotation that depends on mutations
        if col == "ΔΔG":
            X_col = df[ddg_ix, i]
            y_col = labels[ddg_ix]
            if y_col.nunique() < 2:
                results[col] = {'p_value': np.nan, 'log_odds': np.nan, 'std_err': np.nan}
                continue
        else:
            X_col = df[:, i]
            y_col = labels
        
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', ConvergenceWarning)
            X = sm.add_constant(X_col)
            model = sm.Logit(y_col, X)
    
            try:
                result = model.fit(disp=0) 
                p_value = result.pvalues[1] 
                coeff = result.params[1]
                std_err = result.bse[1]
                
            except np.linalg.LinAlgError as e:
                logger.debug("Logistic regression singular matrix: Skipping..")
                logger.debug(e)
                p_value = np.nan
                coeff = np.nan
                std_err = np.nan
                
            except sm.tools.sm_exceptions.PerfectSeparationError as e:
                logger.debug("Logistic regression perfect separation: Skipping..")
                logger.debug(e)
                p_value = np.nan
                coeff = np.nan
                std_err = np.nan

        results[col] = {'p_value': p_value, 'log_odds': coeff, 'std_err': std_err}

    return pd.DataFrame(results)


def uni_log_reg_all_genes(df_annotated, uni_feat_df):
    """
    Univariate logistic regression analysis of all genes.
    """

    results_lst = []
    df_gene_lst = []
    
    for uni_id in df_annotated.Uniprot_ID.unique():
    
        # Process each gene individually
        uni_feat_df_gene = get_uni_feat_for_odds(uni_id=uni_id, uni_feat_df=uni_feat_df)
        
        df_gene = df_annotated[df_annotated["Uniprot_ID"] == uni_id].reset_index(drop=True)
        gene = df_gene.Gene.unique()[0]
        
        df_gene = df_gene.merge(uni_feat_df_gene, on=["Uniprot_ID", "Pos"], how="left")
        target_cols = df_gene.drop(columns=["Pos", "C", "Uniprot_ID", "Gene"]).columns.values
    
        y_data = df_gene["C"]
        X_data = df_gene[target_cols]
        y_data = y_data.fillna(0)
        
        if y_data.nunique() > 1:
            results_gene = uni_log_reg(X_data, y_data)
            results_gene["Gene"] = gene
            results_gene["Uniprot_ID"] = uni_id
            
        else:
            results_gene = pd.DataFrame(np.nan, index=["p_value", "log_odds", "std_err"], columns=target_cols)
            results_gene["Gene"] = gene
            results_gene["Uniprot_ID"] = uni_id
            
        df_gene_lst.append(df_gene)
        results_lst.append(results_gene)

    uni_log_result = pd.concat(results_lst)
    df_genes = pd.concat(df_gene_lst)
    
    return uni_log_result, df_genes


def rename_columns(df):
    """
    Simply rename columns after dummy transformation.
    """
    
    new_columns = {}
    for column in df.columns:
        for prefix in ["SSE", "TYPE"]:
            if column.startswith(f'{prefix}_'):
                new_column = column.replace(f'{prefix}_', '').capitalize()
                new_columns[column] = new_column
                
    return df.rename(columns=new_columns)


def volcano_plot(logreg_results, 
                 fsize=(10, 6), 
                 expand_text_xy=(3, 2), 
                 text_fontsize=10,
                 top_n=15,
                 top_by_gene=False,
                 gene_text=False,
                 save_plot=True,
                 show_plot=False,
                 output_dir=None,
                 cohort="o3d_run"):
    """
    Volcano plot all genes.
    """

    genes = logreg_results[~logreg_results.drop(columns=["Gene", "Uniprot_ID"]).isna().all(axis=1)].Gene.unique()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=DeprecationWarning)
        cmap = plt.cm.get_cmap('tab20', len(genes))
    lgray_rgb = 0.7803921568627451, 0.7803921568627451, 0.7803921568627451, 1.0
    used_colors = []
    
    all_gene_results = []
    plt.figure(figsize=fsize)

    for i, gene in enumerate(genes):
        gene_results = logreg_results[logreg_results["Gene"] == gene].drop(columns=["Gene", "Uniprot_ID"]).dropna(axis=1)
        gene_logodds = gene_results.loc["log_odds", :]
        gene_pvals = gene_results.loc["p_value", :]
        gene_logpvals = -np.log10(gene_pvals)
    
        # Volcano plot
        significant_mask = gene_pvals < 0.01
        non_significant_mask = ~significant_mask
        marker = '*' if cmap(i) in used_colors else 'o'
        used_colors.append(cmap(i))
        plt.scatter(gene_logodds[non_significant_mask], gene_logpvals[non_significant_mask], zorder=1, color='lightgray', alpha=0.7, marker=marker)
        plt.scatter(gene_logodds[significant_mask], gene_logpvals[significant_mask], zorder=2, label=gene, 
                    color="black" if cmap(i) == lgray_rgb else cmap(i), alpha=0.7, marker=marker)
    
        # Append results for annotation
        gene_data = pd.DataFrame({
            'Gene': gene,
            'Log_odds': gene_logodds.values,
            'Pval': gene_pvals.values,
            'Log_pval': gene_logpvals.values,
            'Feature': gene_logodds.index})
        all_gene_results.append(gene_data)
    
    # Annotated top significant n points
    all_gene_results = pd.concat(all_gene_results)
    if top_by_gene:
        top_significant_points = all_gene_results[all_gene_results["Pval"] < 0.01].groupby("Gene").apply(lambda x: x.nsmallest(top_n, 'Pval')).reset_index(drop=True)
    else:
        top_significant_points = all_gene_results[all_gene_results["Pval"] < 0.01].nsmallest(top_n, 'Pval')
    
    annotations = []
    for _, row in top_significant_points.iterrows():
        if gene_text:
            text = f"{row['Gene']}-{row['Feature']}"
        else:
            text = row['Feature']
        annotations.append(plt.text(row['Log_odds'], row['Log_pval'], text, ha='center', va='center', fontsize=text_fontsize, color='black'))
    adjust_text(annotations, expand=expand_text_xy, # expand text bounding boxes by 1.2 fold in x direction and 2 fold in y direction
                arrowprops=dict(arrowstyle='->', color='gray'), lw=0.5)
    
    plt.xlabel('Log odds', fontsize=12)
    plt.ylabel('-log10(p-value)', fontsize=12) 
    plt.axhline(y=-np.log10(0.01), color='lightgrey', linestyle='--', zorder=0)
    plt.axvline(x=0, color='lightgrey', linestyle='--', zorder=0)
    plt.legend(ncol=1 if len(genes) < 20 else 2)
    plt.suptitle(f"{cohort} - Residues' cluster status and annotations associations", y=0.93)    

    if save_plot and output_dir:
        filename = f"{cohort}.volcano_plot.png"
        output_path = os.path.join(output_dir, filename)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.debug(f"Saved {output_path}")
    if show_plot: 
        plt.show()
    plt.close()


def volcano_plot_each_gene(logreg_results,
                           fsize=(3.2, 3),
                           expand_text_xy=(3, 2),
                           text_fontsize=10,
                           top_n=5,
                           ncols=5,
                           all_significant=True,
                           save_plot=True,
                           show_plot=False,
                           output_dir=None,
                           cohort="o3d_run"):
    """
    Volcano plot of individual genes.
    """
            
    genes = logreg_results[~logreg_results.drop(columns=["Gene", "Uniprot_ID"]).isna().all(axis=1)].Gene.unique()
    num_genes = len(genes)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=DeprecationWarning)
        cmap = plt.cm.get_cmap('tab20', num_genes)
    lgray_rgb = 0.7803921568627451, 0.7803921568627451, 0.7803921568627451, 1.0
    
    all_gene_results = []
    
    # Figsize
    if num_genes <= 5:
        ncols = num_genes 
    elif num_genes <= 10:
        ncols = int(np.ceil(num_genes / 2))
    elif num_genes <= 20:
        ncols = 5
    else:
        ncols = 6
    nrows = int(np.ceil(num_genes / ncols))
    
    fsize_x, fsize_y = fsize
    fsize_x = fsize_x * ncols
    fsize_y = fsize_y * nrows
    
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fsize_x, fsize_y), constrained_layout=True)
    # Ensure that the axes is always a 1D subscriptable array
    if not isinstance(axes, Axes):
        axes = axes.flatten()
    axes = np.atleast_1d(axes)
    
    for i, ax in enumerate(axes):
        
        if i < len(genes):
            gene = genes[i]
            gene_results = logreg_results[logreg_results["Gene"] == gene].drop(columns=["Gene", "Uniprot_ID"]).dropna(axis=1)
            gene_logodds = gene_results.loc["log_odds", :]
            gene_pvals = gene_results.loc["p_value", :]
            gene_logpvals = -np.log10(gene_pvals)
        
            # Volcano plot
            significant_mask = gene_pvals < 0.01
            non_significant_mask = ~significant_mask
            ax.scatter(gene_logodds[non_significant_mask], gene_logpvals[non_significant_mask], zorder=1, color='lightgray', alpha=0.7)
            ax.scatter(gene_logodds[significant_mask], gene_logpvals[significant_mask], zorder=2, 
                        color="black" if cmap(i) == lgray_rgb else cmap(i), alpha=0.7)
        
            # Annotated top significant n points
            gene_data = pd.DataFrame({
                'Gene': gene,
                'Log_odds': gene_logodds.values,
                'Pval': gene_pvals.values,
                'Log_pval': gene_logpvals.values,
                'Feature': gene_logodds.index})
            all_gene_results.append(gene_data)
        
            if all_significant:
                top_significant_points = gene_data[gene_data["Pval"] < 0.01]
            else:
                top_significant_points = gene_data.nsmallest(top_n, 'Pval')
            annotations = []
            for _, row in top_significant_points.iterrows():
                annotations.append(ax.text(row['Log_odds'], row['Log_pval'], row['Feature'], ha='center', va='center', fontsize=text_fontsize, color='black'))
            adjust_text(annotations, expand=expand_text_xy, 
                        arrowprops=dict(arrowstyle='->', color='gray'), lw=0.5, ax=ax)
        
            ax.set_xlabel(None)
            ax.set_ylabel(None)
            ax.set_title(gene)
            ax.axhline(y=-np.log10(0.01), color='lightgrey', linestyle='--', zorder=0)
            ax.axvline(x=0, color='lightgrey', linestyle='--', zorder=0)
        else:
            ax.remove()
    
    fig.supxlabel('Log odds')
    fig.supylabel('-log10(p-value)')
    plt.suptitle(f"{cohort} - Residues' cluster status and annotations associations")

    if save_plot and output_dir:
        filename = f"{cohort}.volcano_plot_gene.png"
        output_path = os.path.join(output_dir, filename)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.debug(f"Saved {output_path}")
    if show_plot: 
        plt.show()
    plt.close()
    

def log_odds_plot(logreg_results, 
                  fsize=(1.7,3.6),
                  save_plot=True,
                  show_plot=False,
                  output_dir=None,
                  cohort="o3d_run"):
    """
    Log odds plot.
    """

    genes = logreg_results[~logreg_results.drop(columns=["Gene", "Uniprot_ID"]).isna().all(axis=1)].Gene.unique()
    num_genes = len(genes)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=DeprecationWarning)
        cmap = plt.cm.get_cmap('tab20', num_genes)
    lgray_rgb = 0.7803921568627451, 0.7803921568627451, 0.7803921568627451, 1.0
    
    # Figsize
    if num_genes > 10:
        nrows = 2 
    else:
        nrows = 1
    
    ncols = int(np.ceil(num_genes / nrows))
    fsize_x, fsize_y = fsize
    fsize_x = fsize_x * ncols
    fsize_y = fsize_y * nrows
    
    fig, axes = plt.subplots(nrows, ncols, 
                             figsize=(fsize_x, fsize_y), 
                             sharey=True, 
                             gridspec_kw={'hspace': 0.1*nrows})
    
    # Ensure that the axes is always a 1D subscriptable array
    if not isinstance(axes, Axes):
        axes = axes.flatten()
    axes = np.atleast_1d(axes)

    for i, ax in enumerate(axes):
        
        if i < len(genes):
            gene = genes[i]
            gene_results = logreg_results[logreg_results["Gene"] == gene].drop(columns=["Gene", "Uniprot_ID"])
            gene_logodds = gene_results.loc["log_odds", :]
            gene_pvals = gene_results.loc["p_value", :]
            gene_stderr = gene_results.loc["std_err", :]
        
            # Get 95% confidence interval
            z = 1.96
            lower_ci = np.array(gene_logodds) - z * np.array(gene_stderr)
            upper_ci = np.array(gene_logodds) + z * np.array(gene_stderr)
            
            # Calculate error bars
            lower_error = np.array(gene_logodds) - lower_ci
            upper_error = upper_ci - np.array(gene_logodds)
            error = [lower_error, upper_error]
            
            # Plot
            significant_mask = gene_pvals < 0.01
            non_significant_mask = ~significant_mask
            
            ax.errorbar(gene_logodds, gene_logodds.index.values, yerr=None, xerr=error, fmt='o', capsize=5, capthick=1, markersize=5)
            ax.errorbar(gene_logodds[non_significant_mask], gene_logodds.index.values[non_significant_mask], yerr=None, 
                            xerr=[err[non_significant_mask] for err in error], fmt='o', capsize=5, capthick=1, markersize=5, color='lightgray')
            ax.errorbar(gene_logodds[significant_mask], gene_logodds.index.values[significant_mask], yerr=None, 
                            xerr=[err[significant_mask] for err in error], fmt='o', capsize=5, capthick=1, markersize=5, 
                            color="black" if cmap(i) == lgray_rgb else cmap(i))
            ax.axvline(x=0, color='lightgrey', linestyle='--', zorder=0, lw=1)
            ax.set_xlim(gene_logodds.min()-1.5, gene_logodds.max()+1.5)
            ax.set_xlabel(f"\n\n{gene}", fontsize=12, rotation=0, va='center')
            ax.xaxis.set_label_coords(0.5, 1.11)
        else:
            ax.remove()
        
    plt.suptitle(f"{cohort} - Residues' cluster status and annotations associations", y=1)
    if nrows == 1:
        fig.supxlabel('Log odds', y=-0.015)
        plt.subplots_adjust(top=0.868)
    else:
        fig.supxlabel('Log odds', y=0.037)
        plt.subplots_adjust(top=0.925)


    if save_plot and output_dir:
        filename = f"{cohort}.logodds_plot.png"
        output_path = os.path.join(output_dir, filename)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.debug(f"Saved {output_path}")
    if show_plot: 
        plt.show()
    plt.close()
    
        
def associations_plots(df_annotated, 
                       uni_feat_processed, 
                       output_dir,
                       plot_pars,
                       miss_prob_dict,
                       cohort="o3d_run"):
    """
    Generate volcano plots and log odds plot to look for associations 
    between cluster status and annotated features.
    """
    
    # Prepare data
    df_annotated = df_annotated.copy()
    df_annotated["Miss_prob"] = df_annotated.apply(lambda x: (miss_prob_dict[f"{x.Uniprot_ID}-F{x.F}"][x.Pos-1]), axis=1)
    df_annotated = df_annotated[df_annotated["Miss_prob"] > 0]
    cols_drop = "Mut_in_res", "Mut_in_vol", "Score_obs_sim", "C_ext", "pval", "Clump", "Res", "F", "Ens_Gene_ID", "Ens_Transcr_ID", "PAE_vol", "Miss_prob"
    df_annotated = df_annotated.drop(columns=[col for col in cols_drop if col in df_annotated.columns])
    sse_dummies = pd.get_dummies(df_annotated['SSE'], prefix='SSE')
    df_annotated = pd.concat((df_annotated.drop(columns="SSE"), sse_dummies), axis=1)
    df_annotated = df_annotated.rename(columns={"pLDDT_res" : "pLDDT", "PAE_vol" : "PAE", "DDG" : "ΔΔG"})

    # Add Uniprot features
    uni_feat_processed = uni_feat_processed[uni_feat_processed["Type"] != "REGION"].reset_index(drop=True)
    uni_feat_processed_expanded = expand_uniprot_feat_rows(uni_feat_processed)

    # Perform univariate log reg
    logreg_results, df_annotated_uni_feat = uni_log_reg_all_genes(df_annotated, uni_feat_processed_expanded)
    logreg_results = rename_columns(logreg_results)

    # Plots
    genes = logreg_results[~logreg_results.drop(columns=["Gene", "Uniprot_ID"]).isna().all(axis=1)].Gene.unique()
    if len(genes) > 0:
        output_dir_associations_plots = os.path.join(output_dir, f"{cohort}.associations_plots")
        logger.info(f"Generating associations plots in {output_dir_associations_plots}")
        os.makedirs(output_dir_associations_plots, exist_ok=True)
        log_odds_plot(logreg_results, 
                      output_dir=output_dir_associations_plots, 
                      cohort=cohort,
                      fsize=(plot_pars["log_odds_fsize_x"], plot_pars["log_odds_fsize_y"]))
        volcano_plot(logreg_results, 
                     top_n=plot_pars["volcano_top_n"], 
                     output_dir=output_dir_associations_plots, 
                     cohort=cohort,
                     fsize=(plot_pars["volcano_fsize_x"], plot_pars["volcano_fsize_y"]))
        volcano_plot_each_gene(logreg_results, 
                               output_dir=output_dir_associations_plots, 
                               cohort=cohort,
                               fsize=(plot_pars["volcano_subplots_fsize_x"], plot_pars["volcano_subplots_fsize_y"]))
    else:
        logger.debug("There aren't any relationship to plot: Skipping associations plots..")
    
    uni_feat_cols = ["Uniprot_ID", "Pos", "Domain", "Ptm", "Membrane", "Motif", "Site"]
    df_annotated_uni_feat = df_annotated_uni_feat[[col for col in uni_feat_cols if col in df_annotated_uni_feat.columns]]
    
    return df_annotated_uni_feat


# PLOT WRAPPER
# ============

def generate_plots(gene_result_path,
                   pos_result_path,
                   maf_path,
                   miss_prob_path,
                   seq_df_path,
                   cohort,
                   datasets_dir, 
                   annotations_dir,
                   output_dir,
                   plot_pars,
                   maf_path_for_nonmiss=None,
                   c_genes_only=True,
                   n_genes=30, 
                   lst_genes=None,
                   save_plot=True,
                   show_plot=False,
                   save_csv=True,
                   include_all_pos=False,
                   c_ext=True,
                   title=None,
                   plot_associations=True):
    
    # Load data tracks
    # ================
    
    # Load data
    logger.debug("Loading data")

    gene_result = pd.read_csv(gene_result_path)
    pos_result = pd.read_csv(pos_result_path)
    maf = pd.read_csv(maf_path, sep="\t")
    miss_prob_dict = json.load(open(miss_prob_path))  
    seq_df = pd.read_csv(seq_df_path, sep="\t")    
    uniprot_feat = pd.read_csv(os.path.join(annotations_dir, "uniprot_feat.tsv"), sep="\t")    
    pdb_tool = pd.read_csv(os.path.join(annotations_dir, "pdb_tool_df.tsv"), sep="\t")
    disorder = pd.read_csv(os.path.join(datasets_dir, "confidence.tsv"), sep="\t", low_memory=False)

    # Clean up MOTIF description         TO DO: it should be moved in the build-annotations step
    uniprot_feat.loc[(uniprot_feat["Type"] == "MOTIF") & (
    uniprot_feat["Description"] == "Zinc finger"), "Full_description"] = "Zinc finger"
    uniprot_feat.loc[(uniprot_feat["Type"] == "MOTIF") & (uniprot_feat["Full_description"].str.contains('WIN', case=False)), "Full_description"] = "WIN"
    uniprot_feat.loc[uniprot_feat["Type"] == "MOTIF", "Full_description"] = uniprot_feat.loc[uniprot_feat["Type"] == "MOTIF", "Full_description"].apply(
        lambda x: x.split(";")[0] if len(x.split(";")) > 1 else x)

    # Filter Oncodrive3D result
    logger.debug("Filtering result")
    maf = filter_non_processed_mut(maf, pos_result)
    gene_result, pos_result, genes, uni_ids = filter_o3d_result(gene_result, 
                                                                pos_result, 
                                                                n_genes, 
                                                                lst_genes)
    
    if len(gene_result) > 0:   
        
        # Subset dfs by selected genes and IDs
        logger.debug("Subset genes")
        seq_df, disorder, pdb_tool, uniprot_feat = subset_genes_and_ids(genes, 
                                                                        uni_ids, 
                                                                        seq_df, 
                                                                        disorder, 
                                                                        pdb_tool, 
                                                                        uniprot_feat)

        # Summary plot
        os.makedirs(output_dir, exist_ok=True)
        logger.info(f"Generating summary plot in {output_dir}")
        count_mut_gene_df, count_pos_df, cluster_df = get_summary_counts(gene_result, pos_result, seq_df)
        summary_plot(gene_result, 
                     pos_result, 
                     count_mut_gene_df, 
                     count_pos_df, 
                     cluster_df,
                     output_dir,
                     cohort,
                     plot_pars,
                     save_plot=save_plot,
                     show_plot=show_plot,
                     title=title) 
        
        # Individual gene plots
        if "nonmiss_count" in plot_pars["h_ratios"]:
            maf_nonmiss = get_nonmiss_mut(maf_path_for_nonmiss)
        else:
            maf_nonmiss = None
        output_dir_genes_plots = os.path.join(output_dir, f"{cohort}.genes_plots")
        os.makedirs(output_dir_genes_plots, exist_ok=True)
        logger.info(f"Generating genes plots in {output_dir_genes_plots}")
        
        if c_genes_only:
            n_genes = len(gene_result[gene_result["C_gene"] == 1])
            gene_result, pos_result, genes, uni_ids = filter_o3d_result(gene_result, 
                                                                        pos_result, 
                                                                        n_genes, 
                                                                        lst_genes)
        
        if c_genes_only == False or (c_genes_only and n_genes > 1):
        
            pos_result_annotated, uni_feat_processed = genes_plots(gene_result, 
                                                                    pos_result, 
                                                                    seq_df,
                                                                    maf,
                                                                    maf_nonmiss,
                                                                    miss_prob_dict,
                                                                    output_dir_genes_plots,
                                                                    cohort,
                                                                    annotations_dir,
                                                                    disorder,
                                                                    uniprot_feat,
                                                                    pdb_tool,
                                                                    plot_pars,
                                                                    save_plot=save_plot,
                                                                    show_plot=show_plot,
                                                                    c_ext=c_ext,
                                                                    title=title)
            
            # Associations plots     
            if plot_associations and len(pos_result_annotated) > 0:        
                pos_result_annotated_uni_feat = associations_plots(pos_result_annotated, 
                                                                uni_feat_processed, 
                                                                output_dir,
                                                                plot_pars,
                                                                miss_prob_dict,
                                                                cohort)
                pos_result_annotated = pos_result_annotated.merge(
                    pos_result_annotated_uni_feat.fillna(0), how="left", on=["Uniprot_ID", "Pos"])
                
            
            # Save annotations
            if save_csv and pos_result_annotated is not None:
                logger.info(f"Saving annotated Oncodrive3D result in {output_dir}")
                save_annotated_result(pos_result, 
                                    pos_result_annotated, 
                                    uni_feat_processed, 
                                    output_dir, 
                                    cohort, 
                                    include_all_pos)
                
            logger.info("Plotting completed!")
        else:
            logger.warning("There aren't any significant genes to plot!")
    else:
        logger.warning("There aren't any genes to plot!")
        
        
def generate_comparative_plots(o3d_result_dir_1,
                               cohort_1,
                               o3d_result_dir_2,
                               cohort_2,
                               datasets_dir, 
                               annotations_dir,
                               output_dir,
                               plot_pars,
                               maf_path_nonmiss_1=None,
                               maf_path_nonmiss_2=None,
                               n_genes=30, 
                               lst_genes=None):
    
    # Load data tracks
    # ================
    
    # Load data
    logger.debug("Loading data")
    gene_result_1, pos_result_1, maf_1, miss_prob_dict_1 = load_o3d_result(o3d_result_dir_1, cohort_1)
    gene_result_2, pos_result_2, maf_2, miss_prob_dict_2 = load_o3d_result(o3d_result_dir_2, cohort_2)
    
    seq_df_path = os.path.join(datasets_dir, "seq_for_mut_prob.tsv") 
    seq_df = pd.read_csv(seq_df_path, sep="\t")    
    uniprot_feat = pd.read_csv(os.path.join(annotations_dir, "uniprot_feat.tsv"), sep="\t")    
    pdb_tool = pd.read_csv(os.path.join(annotations_dir, "pdb_tool_df.tsv"), sep="\t")
    disorder = pd.read_csv(os.path.join(datasets_dir, "confidence.tsv"), sep="\t", low_memory=False)
    
    # Clean up MOTIF description         TO DO: it should be moved in the build-annotations step
    uniprot_feat.loc[(uniprot_feat["Type"] == "MOTIF") & (
    uniprot_feat["Description"] == "Zinc finger"), "Full_description"] = "Zinc finger"
    uniprot_feat.loc[(uniprot_feat["Type"] == "MOTIF") & (uniprot_feat["Full_description"].str.contains('WIN', case=False)), "Full_description"] = "WIN"
    uniprot_feat.loc[uniprot_feat["Type"] == "MOTIF", "Full_description"] = uniprot_feat.loc[uniprot_feat["Type"] == "MOTIF", "Full_description"].apply(
        lambda x: x.split(";")[0] if len(x.split(";")) > 1 else x)

    # Filter Oncodrive3D result
    logger.debug("Filtering result")
    maf_1 = filter_non_processed_mut(maf_1, pos_result_1)
    maf_2 = filter_non_processed_mut(maf_2, pos_result_2)
    gene_result_1, pos_result_1, genes_1, uni_ids_1 = filter_o3d_result(gene_result_1, 
                                                                        pos_result_1, 
                                                                        n_genes, 
                                                                        lst_genes)
    gene_result_2, pos_result_2, genes_2, uni_ids_2 = filter_o3d_result(gene_result_2, 
                                                                        pos_result_2, 
                                                                        n_genes, 
                                                                        lst_genes)
    
    # Get shared genes and Uniprot IDs
    shared_genes = [gene for gene in genes_1 if gene in genes_2]
    shared_uni_ids = [uid for uid in uni_ids_1 if uid in uni_ids_2]

    if len(shared_genes) > 0:   
        
        # Subset dfs by selected genes and IDs
        logger.debug("Subset genes")
        seq_df, disorder, pdb_tool, uniprot_feat = subset_genes_and_ids(shared_genes, 
                                                                        shared_uni_ids, 
                                                                        seq_df, 
                                                                        disorder, 
                                                                        pdb_tool, 
                                                                        uniprot_feat)
        
        # Comparative plots
        if "nonmiss_count" in plot_pars["h_ratios"]:
            maf_nonmiss_1 = get_nonmiss_mut(maf_path_nonmiss_1)
            maf_nonmiss_2 = get_nonmiss_mut(maf_path_nonmiss_2)
        else:
            maf_nonmiss_1 = None
            maf_nonmiss_2 = None
        output_dir = os.path.join(output_dir, f"{cohort_1}.{cohort_2}.comparative_plots")
        os.makedirs(output_dir, exist_ok=True)
        
        logger.info(f"Generating comparative plots in {output_dir}")
        comparative_plots(shared_genes,
                        pos_result_1, 
                        maf_1,
                        maf_nonmiss_1,
                        miss_prob_dict_1,
                        cohort_1,
                        pos_result_2,
                        maf_2,
                        maf_nonmiss_2,
                        miss_prob_dict_2,
                        cohort_2,
                        seq_df,
                        output_dir,
                        annotations_dir,
                        disorder,
                        uniprot_feat,
                        pdb_tool,
                        plot_pars,
                        save_plot=True,
                        show_plot=False)
        logger.info("Plotting completed!")
    
    else:
        logger.warning("There aren't any genes to plot!")