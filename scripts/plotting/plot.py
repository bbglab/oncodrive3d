import pandas as pd

import daiquiri
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import os
import json
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
        print(f"Computing missense mut probabilities using mutabilities...")
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
        mut_profile = json.load(open(mutability_config_path))
        print(f"Computing missense mut probabilities...")
        if not isinstance(mut_profile, dict):
            mut_profile = mut_rate_vec_to_dict(mut_profile)
        miss_prob_dict = get_miss_mut_prob_dict(mut_rate_dict=mut_profile, 
                                                seq_df=seq_df)
    else:
        print(f"Mutation profile not provided: Uniform distribution will be used for scoring and simulations.")
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
    directory_path = os.path.join(output_dir, run_name)
    filename = f"{run_name}.summary_plot.png"
    output_path = os.path.join(directory_path, filename)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    logger.info(f"Saved in {output_path}")
    plt.show()
    
      
def genes_plots(gene_result, 
                pos_result, 
                seq_df,
                maf,
                miss_prob_dict,
                count_mut_gene_df, 
                count_pos_df, 
                cluster_df,
                output_dir,
                run_name):   
        

    for gene in gene_result["Gene"].values:

        # Load and parse
        # ==============
        
        uni_id = seq_df[seq_df["Gene"] == gene].Uniprot_ID.values[0]
        af_f = seq_df[seq_df["Gene"] == gene].F.values[0]
        maf_gene = maf[maf["Gene"] == gene]
        mut_count = maf_gene.value_counts("Pos").reset_index()
        pos_result_gene = pos_result[pos_result["Gene"] == gene].sort_values("Pos").reset_index(drop=True)
        pos_result_gene = pos_result_gene[["Pos", "Mut_in_res", "Mut_in_vol", "Ratio_obs_sim", "C", "C_ext", "pval"]]
        pos_result_gene["C"] = pos_result_gene.apply(
            lambda x: 1 if (x["C"] == 1) & (x["C_ext"] == 0) else 2 if (x["C"] == 1) & (x["C_ext"] == 1) else 0, axis=1)
        
        pos_result_gene_hit = pos_result_gene[pos_result_gene["C"] == 1]
        pos_result_gene_not = pos_result_gene[pos_result_gene["C"] == 0]
        pos_result_gene_ext = pos_result_gene[pos_result_gene["C"] == 2]
        pos_vec = pos_result_gene["Pos"].values
        mut_vec = pos_result_gene["Mut_in_res"].values
        vol_vec = pos_result_gene["Mut_in_vol"].values
        max_mut = np.max(pos_result_gene["Mut_in_res"].values)
        max_vol = np.max(pos_result_gene["Mut_in_vol"].values)

        pos_hit = pos_result_gene_hit["Pos"].values
        pos_not = pos_result_gene_not["Pos"].values
        pos_ext = pos_result_gene_ext["Pos"].values
        pos_hit_mut = pos_result_gene_hit["Mut_in_res"].values
        pos_not_mut = pos_result_gene_not["Mut_in_res"].values
        pos_ext_mut = pos_result_gene_ext["Mut_in_res"].values
        pos_hit_vol = pos_result_gene_hit["Mut_in_vol"].values
        pos_not_vol = pos_result_gene_not["Mut_in_vol"].values
        pos_ext_vol = pos_result_gene_ext["Mut_in_vol"].values
        pos_hit_score = pos_result_gene_hit["Ratio_obs_sim"].values
        pos_not_score = pos_result_gene_not["Ratio_obs_sim"].values
        pos_ext_score = pos_result_gene_ext["Ratio_obs_sim"].values
        
        # Get prob vec
        prob_vec = miss_prob_dict[f"{uni_id}-F{af_f}"]   # If none, use uniform
    
        # Get by pos mut, mut_vol, score, etc
        # ===================================
        
        mut_vec = []
        score_vec = []
        #cluster_vec = []
        for pos in range(1, len(prob_vec)+1):

            # Mut count
            if pos in mut_count.Pos.values:
                count = mut_count.loc[mut_count["Pos"] == pos, "count"].values[0]
                score = pos_result_gene.loc[pos_result_gene["Pos"] == pos, "Ratio_obs_sim"].values[0]
            else:
                count = 0
                score = 0
                row_gene = pd.DataFrame({'Pos': [pos], 'Mut_in_res': [0], 'Mut_in_vol': [0], 'Ratio_obs_sim': [np.nan], 'C': [np.nan]})
                pos_result_gene = pd.concat([pos_result_gene, row_gene])
            
            if pos in range(1, len(prob_vec)+1):
                mut_vec.append(count)
                score_vec.append(score)

        pos_result_gene = pos_result_gene.sort_values("Pos").reset_index(drop=True)

        # Normalize score
        mut_freq_vec = mut_vec / sum(mut_vec)
        if np.isnan(score_vec).any():
            score_vec = pd.Series(score_vec).fillna(max(score_vec)).values
        score_norm_vec = score_vec / sum(score_vec) 

        # Generate plot
        # ============= 
        
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, 
                                                 figsize=(24, 12), sharex=True, 
                                                  gridspec_kw={'hspace': 0.05, 
                                                               'height_ratios': [0.22, 0.22, 0.22, 0.32]})
        
        # Plot for Mut_in_res track
        ax1.scatter(pos_hit, pos_hit_mut, label="Significant", color = 'C0', zorder=3)
        ax1.scatter(pos_ext, pos_ext_mut, label="Significant extended", color = 'C2', zorder=2)
        ax1.scatter(pos_not, pos_not_mut, label="Not significant", color = 'C1', zorder=1)
        ax1.fill_between(pos_result_gene['Pos'], 0, max_mut, where=(pos_result_gene['C'] == 1), 
                        color='skyblue', alpha=0.3, label='Mutated *', zorder=0)
        ax1.fill_between(pos_result_gene['Pos'], 0, max_mut, where=((pos_result_gene["C"] == 0) | (pos_result_gene["C"] == 2)), 
                        color='#ffd8b1', alpha=0.6, label='Mutated not *', zorder=0)
        ax1.legend(fontsize=12, ncol=2)

        # Plot for Mut_in_vol track
        ax2.fill_between(pos_result_gene['Pos'], 0, max_vol, where=(pos_result_gene['C'] == 1),
                        color='skyblue', alpha=0.3, label='Mutated *')
        ax2.fill_between(pos_result_gene['Pos'], 0, max_vol, where=((pos_result_gene["C"] == 0) | (pos_result_gene["C"] == 2)), 
                        color='#ffd8b1', alpha=0.6, label='Mutated not *')
        ax2.scatter(pos_not, pos_not_vol, color = 'C1', zorder=1)
        ax2.scatter(pos_hit, pos_hit_vol, color = 'C0', zorder=3)
        ax2.scatter(pos_ext, pos_ext_vol, color = 'C2', zorder=2)

        # Plot for Score track
        ax3.fill_between(pos_result_gene['Pos'], 0, np.max(score_vec), where=(pos_result_gene['C'] == 1), 
                        color='skyblue', alpha=0.3, label='Mutated *')
        ax3.fill_between(pos_result_gene['Pos'], 0, np.max(score_vec), where=((pos_result_gene["C"] == 0) | (pos_result_gene["C"] == 2)), 
                        color='#ffd8b1', alpha=0.6, label='Mutated not *')
        ax3.scatter(pos_hit, pos_hit_score, zorder=3, color="C0")   
        ax3.scatter(pos_not, pos_not_score, zorder=1, color="C1") 
        ax3.scatter(pos_ext, pos_ext_score, zorder=2, color="C2")  

        # Plot for Score and Miss prob track
        max_value = np.max((np.max(mut_freq_vec), np.max(prob_vec), np.max(score_norm_vec)))
        ax4.fill_between(pos_result_gene['Pos'], 0, max_value, where=(pos_result_gene['C'] == 1), 
                        color='skyblue', alpha=0.3, label='Mutated *')
        ax4.fill_between(pos_result_gene['Pos'], 0, max_value, where=((pos_result_gene["C"] == 0) | (pos_result_gene["C"] == 2)), 
                        color='#ffd8b1', alpha=0.6, label='Mutated not *')
        ax4.plot(range(1, len(prob_vec)+1), prob_vec, label="Miss mut prob", zorder=2, color="Red")
        ax4.plot(range(1, len(prob_vec)+1), mut_freq_vec, label="Frequency of mut", zorder=0, color="Purple")
        ax4.plot(range(1, len(prob_vec)+1), score_norm_vec, label="O3D score normalized", zorder=1, color="C2")  
        handles, labels = ax4.get_legend_handles_labels()
        ax4.legend(handles[-3:], labels[-3:], fontsize=12)

        # Set labels and title
        fig.suptitle(f'{gene} - {uni_id}', fontsize=16)
        ax3.set_xlabel('Position', fontsize=14)
        ax1.set_ylabel('Mut in residue', fontsize=14)
        ax2.set_ylabel('Mut in volume', fontsize=14)
        ax3.set_ylabel('O3D score', fontsize=14)
        ax4.set_ylabel('Value', fontsize=14)
        plt.subplots_adjust(top=0.95) 

        # Save
        directory_path = os.path.join(output_dir, run_name)
        filename = f"{run_name}.genes_plot.{gene}_{uni_id}.png"
        output_path = os.path.join(directory_path, filename)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved in {output_path}")
        plt.show()
        
        
# ============
# PLOT WRAPPER
# ============

def generate_plot(gene_result_path, 
                  pos_result_path, 
                  maf, 
                  datasets_dir, 
                  annotations_dir,
                  mut_profile_path,
                  mutability_config_path,
                  output_dir,
                  run_name,
                  n_genes, 
                  significant_only, 
                  lst_genes):
    
    dict_transcripts = {"PTEN" : "ENST00000688308"}              # TO DELETE    
    
    # Load data tracks
    # ================
    
    # Load data
    gene_result = pd.read_csv(gene_result_path)
    pos_result = pd.read_csv(pos_result_path)
    maf = parse_maf_input(maf)
    seq_df_path = os.path.join(datasets_dir, "seq_for_mut_prob.csv") 
    seq_df = pd.read_csv(seq_df_path)    
    pfam = pd.read_csv(os.path.join(annotations_dir, "pfam.tsv"))
    # pfam = pfam.rename(columns={"Ens_transcript_ID" : "Ens_Transcr_ID",        # REMOVE
    #                             "Ens_gene_ID" : "Ens_Gene_ID"})
    pdb_tool = pd.read_csv(os.path.join(annotations_dir, "pdb_tool_df.csv"))
    disorder = pd.read_csv(os.path.join(datasets_dir, "confidence.csv"), low_memory=False)

    # Get processed genes and IDs
    gene_result[gene_result["Status"] == "Processed"].Gene.values
    uni_ids = gene_result.Uniprot_ID.values
    genes = gene_result.Gene.values         
    
    #### HERE IS WHEN I CAN APPLY OTHER FILTERS TO THE GENES & IDs
    
    # Subsett processed genes
    seq_df = seq_df[seq_df["Gene"].isin(genes)]
    seq_df = seq_df[seq_df["Uniprot_ID"].isin(uni_ids)].reset_index(drop=True)
    disorder = disorder[disorder["Uniprot_ID"].isin(uni_ids)].reset_index(drop=True)
    pdb_tool = pdb_tool[pdb_tool["Uniprot_ID"].isin(uni_ids)].reset_index(drop=True)

    # If the Pfam domain is not found for the given transcript, etc.. do this for each key, value in the provided dict
    if dict_transcripts is not None:
        for gene, transcript in dict_transcripts.items():
            seq_df.loc[seq_df["Gene"] == gene, "Ens_Transcr_ID"] = transcript

    # Get subset pfam with gene info
    pfam = seq_df[["Gene", "Uniprot_ID", "Ens_Transcr_ID", "Ens_Gene_ID"]].merge(pfam, how="left", on=["Ens_Transcr_ID", "Ens_Gene_ID"])

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
    
    output_dir_genes_plots = os.path.join(output_dir, "genes_plots")
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
                run_name)
    
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