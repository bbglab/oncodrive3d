import pandas as pd

import daiquiri
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import os
import json
import colorcet as cc
from scripts.run.utils import parse_maf_input

from scripts.run.mutability import init_mutabilities_module
from scripts.run.miss_mut_prob import get_miss_mut_prob_dict, mut_rate_vec_to_dict, get_unif_gene_miss_prob
from scripts import __logger_name__

logger = daiquiri.getLogger(__logger_name__ + ".annotations.plot")



# ==============================
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


# =============
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
        
        
def get_summary_counts(gene_result, pos_result):
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
    count_mut_gene_df = count_mut_gene_df.sort_values("C").reset_index(drop=True)

    # Df with pos count
    pos_result_not = pos_result[pos_result["C"] == 0]
    pos_result_not = pos_result_not.groupby("Gene").apply(len)
    pos_result_not = pos_result_not.reset_index().rename(columns={0 : "Count"})
    pos_result_not["C"] = "Not significant"
    pos_result_hit = pos_result[pos_result["C"] == 1]
    pos_result_hit = pos_result_hit.groupby("Gene").apply(len)
    pos_result_hit = pos_result_hit.reset_index().rename(columns={0 : "Count"})
    pos_result_hit["C"] = "Significant"

    count_pos_df = pd.concat((pos_result_hit, pos_result_not)).sort_values("Gene")
    count_pos_df = count_pos_df.sort_values("C", ascending=False).reset_index(drop=True)

    # Df with cluster count
    cluster_df = pos_result.groupby("Gene").max("Cluster").Cluster.reset_index()
    cluster_df["Cluster"] = cluster_df["Cluster"] + 1
    
    return count_mut_gene_df, count_pos_df, cluster_df
        

def summary_plot(gene_result, 
                 pos_result, 
                 count_mut_gene_df, 
                 count_pos_df, 
                 cluster_df,
                 output_dir,
                 run_name):        

    # Plot
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, 
                                             figsize=(10, 12), 
                                             sharex=True, 
                                             gridspec_kw={'hspace': 0.05, 
                                                          'height_ratios': [0.2, 0.2, 0.2, 0.3]})

    sns.barplot(x='Gene', y='Count', data=count_mut_gene_df, order=gene_result.Gene, hue="C", ax=ax1, ec="black", lw=0.5)
    sns.barplot(x='Gene', y='Count', data=count_pos_df, order=gene_result.Gene, hue="C", ax=ax2,
                palette=sns.color_palette("tab10", n_colors=2), ec="black", lw=0.5)
    sns.barplot(x='Gene', y='Cluster', data=cluster_df, order=gene_result.Gene, ax=ax3, color="lightgray", ec="black", lw=0.5)

    pos_result = pos_result.copy()
    pos_result["C"] = pos_result.C.map({1 : "Significant", 0 : "Not significant"})
    sns.boxplot(x='Gene', y='Ratio_obs_sim', data=pos_result, order=gene_result.Gene, color="lightgray", showfliers=False, ax=ax4)
    sns.stripplot(x='Gene', y='Ratio_obs_sim', data=pos_result, hue="C" ,jitter=True, size=6, alpha=0.5, order=gene_result.Gene.values, 
                palette=sns.color_palette("tab10", n_colors=2), ax=ax4)

    # Add "X" or "O" for significant and not
    for i, gene in enumerate(gene_result['Gene']):
        median_value = pos_result[pos_result['Gene'] == gene]['Ratio_obs_sim'].median()
        is_significant = gene_result[gene_result["Gene"] == gene].C_gene.values[0]
        text_position = median_value - 0.05
        text_symbol = "*" if is_significant else ""
        ax4.text(i, text_position, text_symbol, ha='center', va='center', fontsize=22, fontweight='bold', color="Black")

    # Details
    fig.suptitle(f'Oncodrive3D output', fontsize=16)
    ax1.set_xlabel(None)
    ax2.set_xlabel(None)
    ax4.set_xlabel(None)
    ax1.set_ylabel('Mutations #', fontsize=14)
    ax2.set_ylabel('Mutated residues #', fontsize=14)
    ax3.set_ylabel('Clusters #', fontsize=14)
    ax4.set_ylabel('O3D score\n(obs / μ simulated)', fontsize=14)
    ax1.legend(fontsize=12, loc="upper right")
    ax2.legend(fontsize=12, loc="upper right")
    ax4.legend(fontsize=12)
    plt.xticks(rotation=45, rotation_mode="anchor", ha='right', fontsize=12)
    plt.subplots_adjust(top=0.95) 
    
    # Save
    filename = f"{run_name}.summary_plot.png"
    output_path = os.path.join(output_dir, filename)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    logger.info(f"Saved in {output_path}")
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
    
    maf_nonmiss = pd.read_csv(path_to_maf, sep="\t", dtype={'Chromosome': str})
    maf_nonmiss = maf_nonmiss[maf_nonmiss["Protein_position"] != "-"]
    maf_nonmiss = maf_nonmiss[~(maf_nonmiss['Consequence'].str.contains('Missense_Mutation')
                                | maf_nonmiss['Consequence'].str.contains('missense_variant'))]
    maf_nonmiss = maf_nonmiss[["SYMBOL", 
                               "Consequence", 
                               "Protein_position"]].rename(
                                   columns={"SYMBOL" : "Gene", 
                                            "Protein_position" : "Pos"}).reset_index(drop=True)
                               
    # Parse the consequence with multiple elements
    # Here we assume that the first consequence is the most impactful one
    for i, row in maf_nonmiss.iterrows():
        lst_cnsq = row.Consequence.split(",")
        if len(lst_cnsq) > 1:
            maf_nonmiss.iloc[i].Consequence = lst_cnsq[0]
            
    return maf_nonmiss


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


def check_near_domains(pfam_gene, dist_thr=0.05):
    """
    Check if two domains could be closeR to each other 
    than allowed threshold (ratio of protein size).
    """
        
    near_domains = False
    pfam_gene = pfam_gene.copy()
    pfam_gene = pfam_gene.drop_duplicates(subset='Pfam_name', keep='first')
    mid_pos = (pfam_gene.Pfam_start + pfam_gene.Pfam_end) / 2
    mid_pos_norm = (mid_pos / mid_pos.max()).values
    
    for i in range(len(mid_pos_norm)):
        for j in range(i + 1, len(mid_pos_norm)):
            diff = abs(mid_pos_norm[i] - mid_pos_norm[j])
            if diff < dist_thr:
                near_domains = True

    return near_domains


def genes_plots(gene_result, 
                pos_result, 
                seq_df,
                maf,
                miss_prob_dict,
                count_mut_gene_df, 
                count_pos_df, 
                cluster_df,
                output_dir,
                run_name,
                annotations_dir,
                disorder,
                pfam,
                pdb_tool,
                maf_nonmiss,
                non_missense_count,
                color_cnsq,
                h_ratios,
                s_lw=0.2,
                dist_thr=0.05,
                n_genes=30):   
    """
    Generate a diagnostic plot for each gene showing Oncodrive3D 
    results and annotated features.
    """

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
            
            # Counts
            mut_count = maf_gene.value_counts("Pos").reset_index()
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
            
            # Get prob vec
            prob_vec = miss_prob_dict[f"{uni_id}-F{af_f}"]   # If none, use uniform
        
            # Get by pos mut, score, etc
            # ==========================
            
            score_vec = []
            for pos in range(1, len(prob_vec)+1):

                # Mut count
                if pos in mut_count.Pos.values:
                    score = pos_result_gene.loc[pos_result_gene["Pos"] == pos, "Ratio_obs_sim"].values[0]
                else:
                    score = 0
                    row_gene = pd.DataFrame({'Pos': [pos], 'Mut_in_res': [0], 'Ratio_obs_sim': [np.nan], 'C': [np.nan]})
                    pos_result_gene = pd.concat([pos_result_gene, row_gene])
                
                if pos in range(1, len(prob_vec)+1):
                    score_vec.append(score)

            pos_result_gene = pos_result_gene.sort_values("Pos").reset_index(drop=True)

            # Normalize score
            if np.isnan(score_vec).any():
                score_vec = pd.Series(score_vec).fillna(max(score_vec)).values
            score_norm_vec = np.array(score_vec) / sum(score_vec) 
            
            # Get annotations
            # ===============
            
            disorder_gene = disorder[disorder["Uniprot_ID"] == uni_id].reset_index(drop=True)
            pdb_tool_gene = pdb_tool[pdb_tool["Uniprot_ID"] == uni_id].reset_index(drop=True)
            pfam_gene = pfam[pfam["Uniprot_ID"] == uni_id].reset_index(drop=True)
            ddg = json.load(open(os.path.join(annotations_dir, "stability_change", f"{uni_id}_ddg.json")))
            ddg_vec = avg_per_pos_ddg(pos_result_gene, ddg, maf_gene)
            pos_result_gene["DDG"] = ddg_vec

            # Generate plot
            # ============= 
        
            # Init
            ntracks = len(h_ratios)
            h_ratios_gene = h_ratios.copy()
            ax = 0
            near_domains = check_near_domains(pfam_gene, dist_thr)
            if near_domains: # Make the track for the Pfam domains larger
                h_ratios_gene[len(h_ratios_gene)-1] = 0.1

            # Remove first track if count for non-missense mut is not required
            if non_missense_count == False:
                ntracks -= 1
                del h_ratios_gene[-0]
                h_ratios_gene = np.array(h_ratios_gene) / sum(h_ratios_gene)
                
            fig, axes = plt.subplots(ntracks, 1, 
                                    figsize=(24, 12), sharex=True, 
                                    gridspec_kw={'hspace': 0.1, 
                                                 'height_ratios': h_ratios_gene})

            # Plot for Non-missense mut track   
            # -------------------------------
            if non_missense_count:
                try:
                    if len(mut_count_nonmiss.Consequence.unique()) > 6:
                        ncol = 3
                    else:
                        ncol = 2
                    i = 0
                    axes[ax].vlines(mut_count_nonmiss["Pos"], ymin=0, ymax=mut_count_nonmiss["Count"], color="gray", lw=0.7, zorder=0)
                    for cnsq in mut_count_nonmiss.Consequence.unique():
                        count_cnsq = mut_count_nonmiss[mut_count_nonmiss["Consequence"] == cnsq]
                        if cnsq == "synonymous_variant":
                            order = 1
                        else:
                            order = 2
                        if cnsq in color_cnsq:
                            color = color_cnsq[cnsq]
                        else:
                            color=sns.color_palette("tab10")[i]
                            i+=1
                        axes[ax].scatter(count_cnsq.Pos.values, count_cnsq.Count.values, label=capitalize(cnsq), 
                                        color=color, zorder=order, ec="black", lw=s_lw)   
                    axes[ax].legend(fontsize=11.5, ncol=ncol, framealpha=0.75)
                except Exception as e:
                    logger.warning("Non-missense mutations count not available: Track will be skipped...")
                    logger.warning(f"Encountered the following exception: {e}")
            else:
                ax -= 1

            # Plot for Missense Mut_in_res track
            # ----------------------------------
            axes[ax+1].vlines(mut_count["Pos"], ymin=0, ymax=mut_count["Count"], color="gray", lw=0.7, zorder=1)
            axes[ax+1].scatter(pos_hit, pos_hit_mut, label="Significant", color = 'C0', zorder=4, ec="black", lw=s_lw)   
            axes[ax+1].scatter(pos_ext, pos_ext_mut, label="Significant extended", color = 'C2', zorder=3, ec="black", lw=s_lw)   
            axes[ax+1].scatter(pos_not, pos_not_mut, label="Not significant", color = 'C1', zorder=2, ec="black", lw=s_lw)   
            axes[ax+1].fill_between(pos_result_gene['Pos'], 0, max_mut, where=(pos_result_gene['C'] == 1), 
                            color='skyblue', alpha=0.3, label='Mutated *', zorder=0)
            axes[ax+1].fill_between(pos_result_gene['Pos'], 0, max_mut, where=((pos_result_gene["C"] == 0) | (pos_result_gene["C"] == 2)), 
                            color='#ffd8b1', alpha=0.6, label='Mutated not *', zorder=0)
            axes[ax+1].legend(fontsize=11.5, ncol=2, framealpha=0.75)

            # Plot for Score track
            # --------------------
            axes[ax+2].vlines(pos_result_gene["Pos"], ymin=0, ymax=pos_result_gene["Ratio_obs_sim"], color="gray", lw=0.7, zorder=1)
            axes[ax+2].fill_between(pos_result_gene['Pos'], 0, np.max(score_vec), where=(pos_result_gene['C'] == 1), 
                            color='skyblue', alpha=0.3, label='Mutated *')
            axes[ax+2].fill_between(pos_result_gene['Pos'], 0, np.max(score_vec), where=((pos_result_gene["C"] == 0) | (pos_result_gene["C"] == 2)), 
                            color='#ffd8b1', alpha=0.6, label='Mutated not *')
            axes[ax+2].scatter(pos_hit, pos_hit_score, zorder=3, color="C0", ec="black", lw=s_lw)   
            axes[ax+2].scatter(pos_not, pos_not_score, zorder=1, color="C1", ec="black", lw=s_lw)    
            axes[ax+2].scatter(pos_ext, pos_ext_score, zorder=2, color="C2", ec="black", lw=s_lw)     

            # Plot for Score and Miss prob track
            # ----------------------------------
            max_value = np.max((np.max(prob_vec), np.max(score_norm_vec)))
            axes[ax+3].fill_between(pos_result_gene['Pos'], 0, max_value, where=(pos_result_gene['C'] == 1), 
                            color='skyblue', alpha=0.3, label='Mutated *')
            axes[ax+3].fill_between(pos_result_gene['Pos'], 0, max_value, where=((pos_result_gene["C"] == 0) | (pos_result_gene["C"] == 2)), 
                            color='#ffd8b1', alpha=0.6, label='Mutated not *')
            axes[ax+3].plot(range(1, len(prob_vec)+1), prob_vec, label="Miss mut prob", zorder=2, color="Red")
            axes[ax+3].plot(range(1, len(prob_vec)+1), score_norm_vec, label="O3D score normalized", zorder=1, color="C2")  
            handles, labels = axes[ax+3].get_legend_handles_labels()
            axes[ax+3].legend(handles[-2:], labels[-2:], fontsize=11.5, framealpha=0.75, ncol=2)

            # Plot annotations
            # ================

            # Plot PAE
            # ----------------------------------
            if not np.isnan(pos_result_gene["PAE_vol"]).all():
                max_value = np.max(pos_result_gene["PAE_vol"])
                axes[ax+4].fill_between(pos_result_gene['Pos'], 0, max_value, where=(pos_result_gene['C'] == 1), 
                                color='skyblue', alpha=0.3, label='Mutated *')
                axes[ax+4].fill_between(pos_result_gene['Pos'], 0, max_value, where=((pos_result_gene["C"] == 0) | (pos_result_gene["C"] == 2)), 
                                color='#ffd8b1', alpha=0.6, label='Mutated not *')
                axes[ax+4].plot(pos_result_gene['Pos'], pos_result_gene["PAE_vol"].fillna(0), 
                                label="Confidence", zorder=2, color=sns.color_palette("tab10")[4])

            # Plot disorder
            # -------------
            axes[ax+5].fill_between(pos_result_gene['Pos'], 0, 100, where=(pos_result_gene['C'] == 1), 
                            color='skyblue', alpha=0.3, label='Mutated *')
            axes[ax+5].fill_between(pos_result_gene['Pos'], 0, 100, where=((pos_result_gene["C"] == 0) | (pos_result_gene["C"] == 2)), 
                            color='#ffd8b1', alpha=0.6, label='Mutated not *')
            axes[ax+5].plot(disorder_gene["Pos"], disorder_gene["Confidence"], 
                            label="Confidence", zorder=2, color=sns.color_palette("tab10")[4])

            # Plot pACC
            # ---------
            axes[ax+6].fill_between(pos_result_gene['Pos'], 0, 100, where=(pos_result_gene['C'] == 1), 
                            color='skyblue', alpha=0.3, label='Mutated *')
            axes[ax+6].fill_between(pos_result_gene['Pos'], 0, 100, where=((pos_result_gene["C"] == 0) | (pos_result_gene["C"] == 2)), 
                            color='#ffd8b1', alpha=0.6, label='Mutated not *')
            axes[ax+6].plot(pdb_tool_gene['Pos'], pdb_tool_gene["pACC"].fillna(0), 
                            label="pACC", zorder=2, color=sns.color_palette("tab10")[4])

            # Plot stability change
            # ---------------------
            max_value, min_value = pos_result_gene["DDG"].max(), pos_result_gene["DDG"].min()
            axes[ax+7].fill_between(pos_result_gene['Pos'], min_value, max_value, where=(pos_result_gene['C'] == 1), 
                            color='skyblue', alpha=0.3, label='Mutated *')
            axes[ax+7].fill_between(pos_result_gene['Pos'], min_value, max_value, where=((pos_result_gene["C"] == 0) | (pos_result_gene["C"] == 2)), 
                            color='#ffd8b1', alpha=0.6, label='Mutated not *')
            axes[ax+7].plot(pos_result_gene['Pos'], pos_result_gene["DDG"], 
                            label="Stability change", zorder=2, color=sns.color_palette("tab10")[4])
            axes[ax+7].plot(pos_result_gene['Pos'], np.repeat(0, len(pos_result_gene['Pos'])), 
                            label="No stability change", zorder=2, color='red', linestyle='--', lw=1.2)
            handles, labels = axes[ax+7].get_legend_handles_labels()
            axes[ax+7].legend(handles[-2:], labels[-2:], fontsize=11.5, framealpha=0.75, ncol=2)

            # Clusters label
            # --------------
            clusters_label = pos_result_gene.Cluster.dropna().unique()
            palette = sns.color_palette(cc.glasbey, n_colors=len(clusters_label))
            for i, cluster in enumerate(clusters_label):
                axes[ax+8].fill_between(pos_result_gene['Pos'], -0.5, 0.46, 
                                        where=((pos_result_gene['Cluster'] == cluster) & (pos_result_gene['C'] == 1)),
                                        color=palette[i], alpha=0.6)

            # Secondary structure
            # -------------------
            fill_width = 0.43
            for i, sse in enumerate(['Helix', 'Ladder', 'Coil']):
                c = 0+i
                ya, yb = c-fill_width, c+fill_width
                axes[ax+9].fill_between(pdb_tool_gene["Pos"].values, ya, yb, where=(pdb_tool_gene["SSE"] == sse), 
                                color=sns.color_palette("tab10")[7+i], label=sse)

            # Pfam domains
            # ------------
            pfam_gene = pfam_gene.sort_values("Pfam_start").reset_index(drop=True)
            pfam_color_dict = {}
            
            for n, name in enumerate(pfam_gene["Pfam_name"].unique()):
                pfam_color_dict[name] = f"C{n}"
                
            n = 0
            added_pfam = []
            for i, row in pfam_gene.iterrows():
                if pd.Series([row["Pfam_name"], row["Pfam_start"], row["Pfam_end"]]).isnull().any():
                    continue
                
                name = row["Pfam_name"]
                start = int(row["Pfam_start"])
                end = int(row["Pfam_end"])
                axes[ax+10].fill_between(range(start, end+1), -0.5, 0.48,  
                                alpha=0.5, color=pfam_color_dict[name])
                if name not in added_pfam:
                    if near_domains:
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
                    axes[ax+10].text(((start + end) / 2)+0.5, y, name, ha='center', va='center', fontsize=10, color="black")
                    added_pfam.append(name)

            # Set labels and title
            # --------------------
            fig.suptitle(f'{gene} - {uni_id}', fontsize=16)
            if non_missense_count:
                axes[ax].set_ylabel('Non\nmissense\nmutations', fontsize=13.5)
                axes[ax].set_ylim(-0.5, mut_count_nonmiss["Count"].max()+0.5)

            axes[ax+2].set_xlabel('Position', fontsize=13.5)
            axes[ax+2].set_xlabel(None)
            axes[ax+2].set_ylabel('O3D score', fontsize=13.5)
            axes[ax+1].set_ylabel('Missense\nmutations', fontsize=13.5)
            axes[ax+3].set_ylabel('Value', fontsize=13.5)
            axes[ax+4].set_ylabel('PAE', fontsize=13.5)
            axes[ax+5].set_ylabel('pLDDT', fontsize=13.5)
            axes[ax+5].set_ylim(-10, 110)
            axes[ax+6].set_ylabel('pACC', fontsize=13.5)
            axes[ax+6].set_ylim(-10, 110)
            axes[ax+7].set_ylabel('DDG', fontsize=13.5)
            axes[ax+8].set_ylabel('Clusters             ', fontsize=13.5, rotation=0, va='center')
            axes[ax+8].set_yticks([])  
            axes[ax+8].set_yticklabels([], fontsize=12)
            axes[ax+9].set_yticks([0, 1, 2])  
            axes[ax+9].set_yticklabels(['Helix', 'Ladder', 'Coil'], fontsize=10)
            axes[ax+9].set_ylabel('SSE', fontsize=13.5)
            axes[ax+10].set_yticks([])  
            axes[ax+10].set_yticklabels([], fontsize=12)
            axes[ax+10].set_ylabel('Pfam        ', fontsize=13.5, rotation=0, va='center')
            axes[ax+10].set_ylim(-0.62, 0.6)  

            # Save
            # ----
            filename = f"{run_name}.genes_plot_{j+1}.{gene}_{uni_id}.png"
            output_path = os.path.join(output_dir, filename)
            plt.subplots_adjust(top=0.95) 
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            logger.info(f"Saved in {output_path}")
            plt.close()
        
            
# ============
# PLOT WRAPPER
# ============

def generate_plot(gene_result_path, 
                  pos_result_path, 
                  path_to_maf, 
                  datasets_dir, 
                  annotations_dir,
                  mut_profile_path,
                  mutability_config_path,
                  output_dir,
                  run_name,
                  n_genes=30, 
                  non_significant=False, 
                  lst_genes=None,
                  non_missense_count = False):
    
    dict_transcripts = {"PTEN" : "ENST00000688308"}              # TO DELETE?
    
    
    # XXXXXxxxxx >>>>>>>>>>>>>>>>>>> Parameters <<<<<<<<<<<<<<<< xxxxxXXXXX
    h_ratios = [0.15, 0.15, 0.15, 0.15, 0.1, 0.1, 0.1, 0.1, 0.04, 0.07, 0.04]
    s_lw = 0.2
    dist_thr = 0.05
    
    # Load data tracks
    # ================
    
    # Load data
    gene_result = pd.read_csv(gene_result_path)
    pos_result = pd.read_csv(pos_result_path)
    maf = parse_maf_input(path_to_maf)
    seq_df_path = os.path.join(datasets_dir, "seq_for_mut_prob.csv") 
    seq_df = pd.read_csv(seq_df_path)    
    pfam = pd.read_csv(os.path.join(annotations_dir, "pfam.tsv"))
    pdb_tool = pd.read_csv(os.path.join(annotations_dir, "pdb_tool_df.csv"))
    disorder = pd.read_csv(os.path.join(datasets_dir, "confidence.csv"), low_memory=False)

    # Filter genes and get IDs
    if isinstance(lst_genes, str):
        lst_genes = lst_genes.replace(" ", "")
        lst_genes = lst_genes.split(",")
        gene_result = gene_result[[gene in lst_genes for gene in gene_result["Gene"].values]]    
    if non_significant == False:
        gene_result = gene_result[gene_result["C_gene"] == 1]
    gene_result[gene_result["Status"] == "Processed"].Gene.values
    gene_result = gene_result[:n_genes]                                     
    uni_ids = gene_result.Uniprot_ID.values
    genes = gene_result.Gene.values   
    pos_result = pos_result[[gene in genes for gene in pos_result["Gene"].values]]   
    
    if len(gene_result) > 0:   
    
        # Filter genes in the other df
        seq_df = seq_df[seq_df["Gene"].isin(genes)]
        seq_df = seq_df[seq_df["Uniprot_ID"].isin(uni_ids)].reset_index(drop=True)
        disorder = disorder[disorder["Uniprot_ID"].isin(uni_ids)].reset_index(drop=True)
        pdb_tool = pdb_tool[pdb_tool["Uniprot_ID"].isin(uni_ids)].reset_index(drop=True)

        # If the Pfam domain is not found for the given transcript, etc.. do this for each key, value in the provided dict
        if dict_transcripts is not None:
            for gene, transcript in dict_transcripts.items():
                seq_df.loc[seq_df["Gene"] == gene, "Ens_Transcr_ID"] = transcript

        # Get subset pfam with gene info
        pfam = seq_df[["Gene", "Uniprot_ID", "Ens_Transcr_ID", "Ens_Gene_ID"]].merge(
            pfam, how="left", on=["Ens_Transcr_ID", "Ens_Gene_ID"])

        # Get missense mut probability dict
        miss_prob_dict = get_miss_mut_prob_for_plot(mut_profile_path, mutability_config_path, seq_df)
        
        # Summary plots
        # =============
        
        create_plot_dir(output_dir)
        count_mut_gene_df, count_pos_df, cluster_df = get_summary_counts(gene_result, pos_result)
        summary_plot(gene_result, 
                    pos_result, 
                    count_mut_gene_df, 
                    count_pos_df, 
                    cluster_df,
                    output_dir,
                    run_name) 
        
        # Plots for individual genes
        # ==========================
        
        # Get non missense mutations
        if non_missense_count:
            maf_nonmiss = get_nonmiss_mut(path_to_maf)
        else:
            maf_nonmiss = None
        
        # Cnsq color
        color_cnsq = {"start_lost" : "C1",
                      "start_lost" : "C5",
                      "stop_gained" : "C3",
                      "stop_lost" : "C4",
                      "splice_region_variant" : "C2",
                      "splice_acceptor_variant" : "C6",
                      "splice_donor_variant" : "C7",
                      "stop_retained_variant" : "C8",
                      "synonymous_variant" : "C9"}
        
        output_dir_genes_plots = os.path.join(output_dir, f"{run_name}.genes_plots")
        create_plot_dir(output_dir_genes_plots)
        genes_plots(gene_result, 
                    pos_result, 
                    seq_df,
                    maf,
                    miss_prob_dict,
                    count_mut_gene_df, 
                    count_pos_df, 
                    cluster_df,
                    output_dir_genes_plots,
                    run_name,
                    
                    annotations_dir,
                    disorder,
                    pfam,
                    pdb_tool,
                    maf_nonmiss,
                    non_missense_count,
                    color_cnsq,
                    h_ratios,
                    s_lw,
                    dist_thr,
                    n_genes)
        
        logger.info("Plotting completed!")
        
    else:
        logger.info("There aren't any genes to plot!")
    
    ## Select the genes from gene_result
    # IF: processed only
    # IF: hits only
    # IF: top
    
    ## Subset seq_df
    
    ## Get miss mut prob 
    
    ## Generate directory output + output/gene_plots

    ## Retrieve all annotations for selected genes    
    # Disorder
    # PAE
    # DDG (average across mut)
    # Secondary structure
    # pAAC
    
    ## Summary plot
    ## Summary plot extra (extra if there is time)
    
    ## Diagnostic plots