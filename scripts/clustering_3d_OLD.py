""" 
The module includes the functions and the script to run 
an HotMAPs-inspired method that for each gene uses only 
the canonical predicted structure found in AlphaFold db

###################################### EXAMPLE USAGE ############################################

python3 clustering_3d.py -i /workspace/datasets/intogen/runs/20210108/intogen_20210108/steps/hotmaps/PCAWG_WGS_COLORECT_ADENOCA.in.maf \           
-o /workspace/projects/alphafold_features/hotmaps/test_result \
-P ../required_files/extra/mut_profile/prior.weighted.normalized.PCAWG_WGS_COLORECT_ADENOCA.json \
-s ../required_files/seq_for_mut_prob.csv \
-c ../required_files/cmaps/ \
-u ../required_files/af_uniprot_to_gene_id.json \
-n 10000 \
-t COREAD \
-C PCAWG_WGS_COLORECT_ADENOCA.V1.2


-p 

#################################################################################################
"""


from lib2to3.pytree import generate_matches
import os
import json
import numpy as np
import seaborn as sns
import pandas as pd
import re
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from math import sqrt
from math import floor
import statistics
import pickle
import re
import scipy as sp
from Bio.Seq import Seq
from Bio import SeqIO
from itertools import product
from difflib import SequenceMatcher
from progressbar import progressbar
from tqdm.notebook import tqdm
import argparse
import networkx as nx
import networkx.algorithms.community as nx_comm
from scipy.stats import rankdata
import glob
import time
from utils.utils import get_pos_fragments, uniprot_to_hugo
from utils.miss_mut_prob import get_miss_mut_prob_dict, mut_rate_vec_to_dict



## Load data

def get_uniprot_from_gene(gene_name, dictio):
    """
    Simple function that, given a dictionary of gene name 
    (keys) and Uniprot IDs (values), and a specific gene 
    name, returns the corresponding Uniprot ID
    """
    
    return list({k:v for k, v in dictio.items() if v == gene_name})[0]
    

def parse_maf_input(maf_input_path):
    """
    Parse in.maf file which is used as 
    input for the HotMAPS method.
    """

    # Load
    maf = pd.read_csv(maf_input_path, sep="\t")

    # Select only missense mutation and extract Gene name and mut
    maf = maf.loc[maf.Variant_Classification == "Missense_Mutation"].copy()
    maf["Pos"] = maf.loc[:, "HGVSp_Short"].apply(lambda x: int(re.sub("\\D", "", (x[2:]))))
    maf["WT"] = maf["HGVSp_Short"].apply(lambda x: re.findall("\\D", x[2:])[0])
    maf["Mut"] = maf["HGVSp_Short"].apply(lambda x: re.findall("\\D", x[2:])[1])
    maf = maf[["Hugo_Symbol", "Pos", "WT", "Mut"]]
    maf = maf.sort_values("Pos").rename(columns={"Hugo_Symbol" : "Gene"})
    
    return maf


def parse_cluster_output(out_cluster_path):
    """
    Parse out.gz file which is the output of HotMAPS.
    It will be used as ground truth for the new clustering method.
    """
    
    # Select only necessary info and rename to match the input file
    cluster = pd.read_csv(out_cluster_path, sep="\t")
    cluster = cluster.copy().rename(columns={"CRAVAT Res" : "Pos"})[["Pos", "Ref AA", "HUGO Symbol", "Min p-value", "q-value"]]
    cluster = cluster.rename(columns={"HUGO Symbol" : "Gene", "Ref AA" : "WT"})
    
    return cluster



## Eval


def get_evaluation_df(test_result, mut_gene_df, drop_extra=False):
    """
    Merge the results from the statistical test (excess of densities 
    of mutations) with the score assigned from HotMAPs.
    """

    if drop_extra:
        cols_to_drop = ["Gene", "WT", "Mut"]
        mut_gene_df = mut_gene_df.drop(columns = cols_to_drop)
    
    # Merge the new result (mutations predicted to be in positions enriched of mutations)
    # with the original data (mutations that were not predicted to be enriched will have NA in that col)
    mut_gene_df = mut_gene_df.sort_values("Pos")
    evaluation_df = test_result.merge(mut_gene_df, how = "outer").sort_values("Pos").reset_index(drop=True)
    
    return evaluation_df


def evaluate(evaluation_df):

    evaluation_df = evaluation_df[["Gene", 
                                   "Pos", 
                                   "New_hotmaps", 
                                   "HotMAPs"]].drop_duplicates(keep="first")
    
    all_mut = len(evaluation_df)

    all_pos = evaluation_df[evaluation_df["HotMAPs"] == 1]
    pos_pred = evaluation_df[evaluation_df["New_hotmaps"] == 1]
    tp_pred = pos_pred[pos_pred["HotMAPs"] == 1]
    fp_pred = pos_pred[pos_pred["HotMAPs"] == 0]

    all_neg = evaluation_df[evaluation_df["HotMAPs"] == 0]
    neg_pred = evaluation_df[evaluation_df["New_hotmaps"] == 0]
    tn_pred = neg_pred[neg_pred["HotMAPs"] == 0]
    fn_pred = neg_pred[neg_pred["HotMAPs"] == 1]
    
    if len(all_pos) > 0:
        tpr = round(len(tp_pred) / (len(tp_pred) + len(fn_pred)), 3)
        if len(pos_pred) > 0:
            precision = round(len(tp_pred) / (len(tp_pred) + len(fp_pred)), 3)
        else:
            precision = np.nan
    else:
        tpr = np.nan
        precision = np.nan
    if len(all_neg):
        tnr = round(len(tn_pred) / (len(tn_pred) + len(fp_pred)), 3)
    else:
        tnr = np.nan

    return pd.DataFrame({"Gene" : evaluation_df.Gene.unique()[0],
                         "All_mut" : all_mut, 
                         "All_pos" : len(all_pos),
                         "Pos_pred" : len(pos_pred), 
                         "TP_pred" : len(tp_pred), 
                         "FP_pred" : len(fp_pred), 
                         "All_neg" : len(all_neg),
                         "Neg_pred" : len(neg_pred),
                         "TN_pred" : len(tn_pred),
                         "FN_pred" : len(fn_pred),
                         "TPR" : tpr, 
                         "Precis" : precision,
                         "TNR" : tnr,}, 
                        index=[1])

    
def plot_venn(eval_df_final, save=False, path="venn.png"):
    """
    Venn diagram for genes found by HotMAPS, 
    the New HotMAPS, and both methods.
    """
    
    both = len(eval_df_final[(eval_df_final.New_hotmaps == 1) & (eval_df_final.HotMAPs == 1)].Gene.unique())
    new_hotmaps = len(eval_df_final[eval_df_final.New_hotmaps == 1].Gene.unique()) - both
    hotmaps = len(eval_df_final[eval_df_final.HotMAPs == 1].Gene.unique()) - both

    venn2(subsets = (hotmaps, new_hotmaps, both), set_labels = ('HotMAPS', 'New HotMAPS'))
    plt.title("Genes with one or more clusters", fontsize=14)
    if save: 
        plt.savefig(path, bbox_inches="tight", dpi=300)
        plt.clf()
    else:
        plt.show()


def fdr(p_vals):
    """
    Compute false discovery rate.
    """

    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1

    return fdr


def sort_by_pval_and_mut(result_gene):
    """
    Sort by p-value and break the tie using the 
    density of mutations (larger the better) and 
    total mutation in the gene (smaller the better). 
    """
    
    result_gene = result_gene.copy()
    result_gene["Mut_in_vol"] = -result_gene["Mut_in_vol"]
    result_gene = result_gene.sort_values(["pval", "Mut_in_vol", "Mut_count"], ascending=True).reset_index(drop=True)
    result_gene["Mut_in_vol"] = -result_gene["Mut_in_vol"]
    
    return result_gene


def get_final_gene_result(result_pos, result_gene, alpha_gene=0.05):
    """
    Output the final dataframe including gene global pval, 
    significant positions, communities, processing status, etc.
    """
    
    # Get significant positions and communities for each gene
    clusters = result_pos[result_pos["C"] == 1].groupby("Gene").apply(lambda x: (x["Pos"].values)).reset_index().rename(
        columns={0 : "C_pos"})
    clusters["C_community"] = result_pos[result_pos["C"] == 1].groupby("Gene").apply(lambda x: x["Community"].values).reset_index(drop=True)
    
    # Annotate each gene with significant hits
    result_gene = clusters.merge(result_gene, on="Gene", how="outer")
    
    # Get gene pval, qval, and largest density omong the most significant hits
    gene_pvals = result_pos.groupby("Gene").apply(lambda x: min(x["pval"].values)).reset_index().rename(
        columns={0 : "pval"})
    gene_pvals["qval"] = fdr(gene_pvals["pval"])
    gene_pvals["Mut_in_vol"] = result_pos.groupby("Gene").apply(lambda x: max(x[x["pval"] == min(x["pval"])].Mut_in_vol)).values

    # Combine pval to significant annotation
    result_gene = gene_pvals.merge(result_gene, on="Gene", how="outer")

    # Add gene binary clustering label and sort genes
    gene_label = result_gene.apply(lambda x: 1 if x.qval < alpha_gene else 0, axis=1)
    result_gene.insert(3, "C_gene", gene_label)
    result_gene.insert(7, "Mut_in_vol", result_gene.pop("Mut_in_vol"))
    result_gene = sort_by_pval_and_mut(result_gene)
    
    return result_gene


## Clustering score and simulations


def get_anomaly_score(vec_mut_in_vol, gene_mut, vec_vol_miss_mut_prob):          
    """
    Compute a metric that scores the anomaly of observing a certain 
    number of mutations in the volume of a residue.
    It takes into account the volume and the mutation rate of the codon 
    of each residue within that volume.
    
    Score: loglik equal or larger mut_count / loglik(N)
    """
    
    den = sp.stats.binom.logpmf(k=gene_mut, n=gene_mut, p=vec_vol_miss_mut_prob)
    return sp.stats.binom.logsf(k=vec_mut_in_vol-1, n=gene_mut, p=vec_vol_miss_mut_prob) / den


def simulate_mutations(n_mutations, p, size):
    """
    Simulate the mutations given the mutation rate of a cohort.
    """

    rng = np.random.default_rng()
    samples = rng.multinomial(n_mutations, p, size=size)
    return samples


def get_sim_anomaly_score(mut_count, 
                            cmap,
                            gene_miss_prob,
                            vol_missense_mut_prob,
                            num_iteration=1000):
    """
    Simulated mutations following the mutation profile of the cohort.
    Compute the log-likelihood of observing k or more mutation in the 
    volume and compare it with the corresponding simualted rank.
    """
    
    # Generate x sets of random mutation distributed following the mut 
    # profile of the cohort, each with the same size of the observed mut   
    mut_sim = simulate_mutations(mut_count, gene_miss_prob, num_iteration)
    
    # Get the density of mutations of each position in each iteration
    density_sim = np.einsum('ij,jk->ki', cmap, mut_sim.T.astype(float), optimize=True)
    
    # Compute the ranked score of the densities obtained at each iteration
    # sign is used to sort in descending order
    loglik_plus = -np.sort(-get_anomaly_score(density_sim, mut_count, vol_missense_mut_prob))
    
    return pd.DataFrame(loglik_plus).T



## Communities detection 


def get_network(nodes, mut_count_v, cmap):
    """
    Generate a network with significative pos as nodes 
    and ratio of shared mutation (Jaccard score) as edges. 
    """
    
    G = nx.Graph()
    G.add_nodes_from(nodes)

    # Iterate through each hit
    for i, ipos in enumerate(nodes):
        for j, jpos in enumerate(nodes):
            if i > j:
                
                # Add an edge if they are in contact
                ix, jx = ipos-1, jpos-1
                if cmap[ix, jx] == 1:

                    # Get the common res and their mut
                    neigh_vec_i, neigh_vec_j = cmap[ix], cmap[jx]
                    common_neigh = neigh_vec_j * neigh_vec_i
                    num_mut = np.dot(common_neigh, mut_count_v)

                    # Get the sum of the union of the mut
                    all_neigh = (neigh_vec_i + neigh_vec_j != 0).astype(int)
                    union_num_mut = np.dot(all_neigh, mut_count_v)

                    # Compute the Jaccard score or avg ratio between ij shared mut    
                    jaccard = np.round(num_mut/union_num_mut, 3)
                    G.add_edge(ipos, jpos, weight = jaccard)
    return G


def get_community_index_nx(pos_hits, communities):
    
    """
    Parse the labels returned by communities detection algorithms from NetworkX.
    """
    
    communities_mapper = {}
    for ic, c in enumerate(communities):
        for p in c:
            communities_mapper[p] = ic

    return pos_hits.map(communities_mapper).values



## Wrapper function
def pseudo_hotmaps(gene, 
                   mut_gene_df,                                     
                   gene_to_uniprot_dict, 
                   cmap_path,
                   miss_prob_dict,
                   fragment=1,
                   alpha=0.01,
                   num_iteration=10000):
    """
    Compute local density of missense mutations for a sphere of 10A around          
    each amino acid position of the selected gene product. It performed a 
    rank-based comparison between observed density and simulated ones in 
    absense of positive selection (cohort mut profile). Get an experimental 
    p-val for the local enrichment.
    
    Parameters:                                                                       ##  CHANGE/UPDATE
    -----------
    gene : str
    mut_gene_df : pandas dataframe
        It must include the mutated positions of the gene. 

    gene_to_uniprot_dict : dict
        Uniprot_ID as key and corresponding genes as values
        
    neighbours_df : pandas df
    miss_prob_dict : dict
        Uniprot_ID as keys and per residue prob of missense mut as values
        
    af_structures_path : str
    num_iteration : int
    v : bolean
    plot_contact_map : bolean

    Returns:                                                                            ##  CHANGE/UPDATE
    ----------
    evaluation_df : pandas df (per-position evaluation)
    test_result : pandas df (per-gene summary evaluation)
    status_df : pandas df (per-gene processing status)
    """
    

    ## Initialize
    mut_count = len(mut_gene_df)
    result_gene_df = pd.DataFrame({"Gene" : gene,
                                   "No_frag" : fragment,                           
                                   "Mut_count" : mut_count,
                                   "Max_mut_pos" : np.nan,
                                   "Structure_max_pos" : np.nan,
                                   "Status" : np.nan}, 
                                    index=[1])
    
    # Load cmap
    cmap_complete_path = f"{cmap_path}/{gene_to_uniprot_dict[gene]}-F{fragment}.npy"
    if os.path.isfile(cmap_complete_path):
        cmap = np.load(cmap_complete_path) 
    else:
        result_gene_df["Status"] = "Cmap_not_found"
        return None, result_gene_df

    # Check if there is a mutation that is not in the structure      
    if max(mut_gene_df.Pos) > len(cmap):
        result_gene_df["Max_mut_pos"] = max(mut_gene_df.Pos)  
        result_gene_df["Structure_max_pos"] = len(cmap)
        result_gene_df["Status"] = "Mut_not_in_structure"
        return None, result_gene_df


    ## Get expected local myssense mutation density

    # Probability that each residue can be hit by a missense mut
    gene_miss_prob = np.array(miss_prob_dict[f"{gene_to_uniprot_dict[gene]}-F{fragment}"])

    # Probability that the volume of each residue can be hit by a missense mut
    vol_missense_mut_prob = np.dot(cmap, gene_miss_prob)
    
    
    ## Get observed and ranked simulated scores (loglik+_LFC)
    
    # Get the observed mut count and densities 
    count = mut_gene_df.Pos.value_counts()
    mut_count_v = np.zeros(len(cmap))
    mut_count_v[count.index - 1] = count.values
    mut_count_m = mut_count_v.reshape((1, -1))
    density_m = np.einsum('ij,jk->ki', cmap, mut_count_m.T, optimize=True)
    mutated_pos = np.sort(count.index)

    # Do not process if there isn't any density larger than 1
    if max(density_m[0][mutated_pos-1]) <= 1:                       
        result_gene_df["Status"] = "No_density"
        return None, result_gene_df

    if max(density_m[0][mutated_pos-1]) > 10:                                    ################
        print(gene, max(density_m[0][mutated_pos-1]))
    
    # Inialize result df 
    result_pos_df = pd.DataFrame({"Pos" : mutated_pos, "Mut_in_vol" : density_m[0, mutated_pos-1].astype(int)})

    # Get the ranked simulated score
    sim_anomaly = get_sim_anomaly_score(mut_count, 
                                        cmap, 
                                        gene_miss_prob,
                                        vol_missense_mut_prob,                                                                               
                                        num_iteration=num_iteration) 

    # Get ranked observed score (loglik+_LFC) 
    no_mut_pos = len(result_pos_df)
    sim_anomaly = sim_anomaly.iloc[:no_mut_pos,:].reset_index()
    result_pos_df["Abs_anomaly"] = get_anomaly_score(result_pos_df["Mut_in_vol"], mut_count, vol_missense_mut_prob[result_pos_df["Pos"]-1])
    result_pos_df = result_pos_df.sort_values("Abs_anomaly", ascending=False).reset_index(drop=True)
    result_pos_df["Avg_sim_anomaly"] = sim_anomaly.iloc[:,1:].mean(axis=1)
   # result_pos_df["Diff_obs_avg_sim"] = result_pos_df["Abs_anomaly"] - result_pos_df["Avg_sim_anomaly"]

    ## Compute p-val and assign hits

    # Add to the simulated score of each iteration its standard deviation  
    # This makes the method more conservative (avoid borderline cases)
    sim_anomaly.iloc[:,1:] = sim_anomaly.apply(lambda x: x[1:] + x[1:].std(), axis=1)

    # Experimental p-val
    result_pos_df["pval"] = sim_anomaly.apply(lambda x: sum(x[1:] >= result_pos_df["Abs_anomaly"].values[int(x["index"])]) / len(x[1:]), axis=1)
    # result_pos_df["n_sim_≥_obs"] = sim_anomaly.apply(lambda x: sum(x[1:] >= result_pos_df["Abs_anomaly"].values[int(x["index"])]), axis=1)
    # result_pos_df["n_sim"] = sim_anomaly.apply(lambda x: len(x[1:]), axis=1)

    # Assign hits
    result_pos_df["C"] = [int(i) for i in result_pos_df["pval"] < alpha]                

    
    ## Communities detection
    
    pos_hits = result_pos_df[result_pos_df["C"] == 1].Pos
    if len(pos_hits) > 0:
        if len(pos_hits) > 1:
            # Build network and perform detection
            G = get_network(pos_hits, mut_count_v, cmap)
            communities = nx_comm.label_propagation_communities(G)
            clusters = get_community_index_nx(pos_hits, communities)

        else:
            # Assign cluster 0 to the only pos hit
            clusters = 0 
        meta_clusters = pd.DataFrame({"Pos" : pos_hits, "Community" : clusters})
        result_pos_df = result_pos_df.merge(meta_clusters, how = "left", on = "Pos")
    else:
        result_pos_df["Community"] = np.nan
        

    ## Evaluate the result and output
    result_pos_df.insert(0, "Gene", gene)
    result_pos_df.insert(2, "Mut_in_gene", mut_count)
    result_gene_df["Status"] = "Processed"
    
    return result_pos_df, result_gene_df



#######################

# STEP 1

# MODIFY THE FOLLOWING FUNCTION SO THAT
#   - Run for one cohort   
#   - Takes as input the paths of a maf, a mut_rate, output dir and name
#   - Outputs maf with clustering info, status processing, log_file

#######################

# STEP 2 (only if there is time)

# MODIFY THE STRUCTURE OF THE SCRIPTS, DIVIDE IT INTO MULTIPLE FILES IN DIFFERENT FOLDERS
#   - Keep a main in the 3d_clustering folder (wrapper function)
#   - Name 3d_clustering.py the script that includes the function to run clustering for one gene

#######################

# STEP 3 

# Use a docker to allow anyone to use the easily use the tool

#######################

# STEP 4 (only if there is time)

# Convert big functions into classes



######################### EXAMPLE ##########################

# -m /workspace/datasets/intogen/runs/20210108/intogen_20210108/steps/hotmaps/PCAWG_WGS_COLORECT_ADENOCA.in.maf"
# -p /workspace/projects/alphafold_features/hotmaps/mut_rate/miss_mut_prob/PCAWG_WGS_COLORECT_ADENOCA.miss_mut_prob.json
# -i 100
# -o /workspace/projects/alphafold_features/hotmaps/test_result
# -t CESC
# -c HARTWIG_UTERUS_CERVICAL


# "args": ["-m /workspace/datasets/intogen/runs/20210108/intogen_20210108/steps/hotmaps/PCAWG_WGS_COLORECT_ADENOCA.in.maf",
#          "-p /workspace/projects/alphafold_features/hotmaps/mut_rate/miss_mut_prob/PCAWG_WGS_COLORECT_ADENOCA.miss_mut_prob.json",
#          "-i 100",
#          "-o /workspace/projects/alphafold_features/hotmaps/test_result",
#          "-t COREAD",
#          "-c PCAWG_WGS_COLORECT_ADENOCA"]

############################################################



def main():
    """
    Wrapper function.
    """

    ## Initialize
    version = "2023_v1.2"    # LAST CHANGE: restructured pipeline

    # Parser
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)

    parser.add_argument("-i", "--input_maf", help="Path of the maf file used as input", type=str, required=True)
    parser.add_argument("-o", "--output_dir", help="Path to output directory", type=str, required=True)

    group.add_argument("-p", "--miss_mut_prob", help="Path to the dict of missense mut prob of each protein based on mut profile of the cohort (json)",  type=str)
    group.add_argument("-P", "--mut_profile", help="Path to the mut profile (list of 96 floats) of the cohort (json)", type=str)
    parser.add_argument("-s", "--seq_df", 
                        help="Path to the dataframe including DNA and protein seq of all gene/proteins (all AF predicted ones)", 
                        type=str, 
                        default="../required_files/seq_for_mut_prob.csv")       

    parser.add_argument("-c", "--cmap_path", help="Path to the directory containting the contact map of each protein", type=str, default="../required_files/cmaps/")
    parser.add_argument("-u", "--uniprot_to_gene_dict", help="Path to a dictionary including Uniprot_ID : HUGO symbol mapping", type=str)  

    parser.add_argument("-n", "--n_iterations", help="Number of densities to be simulated", type=int, default=10000)
    parser.add_argument("-a", "--alpha_level_res", help="Significant threshold for the p-value of protein residues", type=float, default=0.01)
    parser.add_argument("-A", "--alpha_level_gene", help="Significant threshold for the global q-value of the genes", type=float, default=0.05)
    
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
    uniprot_to_gene_path = args.uniprot_to_gene_dict

    print(f"Starting 3D-clustering [{version}]..\n")
    print(f"Iterations: {num_iteration}")
    print(f"Significant level (position): {alpha_pos}")
    print(f"Global significant level (gene): {alpha_gene}")
    print(f"Cohorts: {cohort}")
    print(f"Cancer type: {cancer_type}")
    print(f"Output directory: {output_dir}")


    ## Load necessary files

    # Uniprot to HUGO symbol mapping
    if uniprot_to_gene_dict is not None:
        uniprot_to_gene_dict = json.load(open(uniprot_to_gene_path))
        gene_to_uniprot_dict = {g:p for (p,g) in uniprot_to_gene_dict.items()}
    else:
        uniprot_to_gene_dict = uniprot_to_hugo()

    # MAF input
    data = parse_maf_input(maf_input_path)

    # Missense mut prob  
    if miss_mut_prob_path is not None:
        # Load dict with miss prob of each prot
        miss_prob_dict = json.load(open(miss_mut_prob_path))
    else:
        # Compute dict from mut profile of the cohort and dna sequences
        mut_profile = json.load(open(mut_profile_path))
        seq_df = pd.read_csv(seq_df_path)
        print("\nComputing missense mut probability for each protein")
        mut_profile = mut_rate_vec_to_dict(mut_profile)
        miss_prob_dict = get_miss_mut_prob_dict(mut_rate_dict=mut_profile, seq_df=seq_df)


    ## Run

    result_pos_lst = []
    result_gene_lst = []

    # Get genes with enough mut
    genes = data.groupby("Gene").apply(lambda x: len(x))
    genes_mut = genes[genes >= 2].index
    genes_no_mut = genes[genes < 2].index
    result_gene = pd.DataFrame({"Gene" : genes_no_mut,
                            "No_frag" : np.nan,
                            "Mut_count" : 1,
                            "Max_mut_pos" : np.nan,
                            "Structure_max_pos" : np.nan,
                            "Status" : "No_mut"})
    result_gene_lst.append(result_gene)

    # Process gene
    print("Performing 3D clustering on genes with enough mutations..")
    for gene in progressbar(genes_mut):

        ########  XXXXX ########

        mut_gene_df = data[data["Gene"] == gene]

        # If there is a single fragment
        if seq_df[seq_df["Gene"] == gene].F.max() == 1:      

            # Run clustering
            try:
                pos_result, result_gene = pseudo_hotmaps(gene, 
                                                        data,
                                                        gene_to_uniprot_dict, 
                                                        cmap_path,
                                                        miss_prob_dict,
                                                        fragment=1,
                                                        alpha=alpha_pos,
                                                        num_iteration=num_iteration)
                result_gene_lst.append(result_gene)
                
                if pos_result is not None:
                    result_pos_lst.append(pos_result)
        
            # Can't run clustering                                                                >>>> Should raise a better exception to capture a more specific error
            except:
                result_gene = pd.DataFrame({"Gene" : gene,
                                        "No_frag" : np.nan,
                                        "Mut_count" : np.nan,
                                        "Max_mut_pos" : np.nan,
                                        "Structure_max_pos" : np.nan,
                                        "Status" : "Not_processed"},
                                        index=[0])
                result_gene_lst.append(result_gene)

        # If the protein is fragmented
        else:
            mut_gene_df["F"] = get_pos_fragments(mut_gene_df)
            f_pos_result_lst = []
            f_result_gene_lst = []
            prot_community = 0

            # Consider each fragment as a individual protein
            for fragment in mut_gene_df["F"].unique():

                mut_fragment_df = mut_gene_df[mut_gene_df["F"] == fragment]

                # Use relative pos of the fragments
                mut_fragment_df = mut_fragment_df.copy()
                mut_fragment_df["Pos"] = mut_fragment_df["Pos"] - 1400 * (fragment-1)

                # Perform clustrering on fragment
                f_pos_result, f_result_gene = pseudo_hotmaps(gene, 
                                                         mut_fragment_df,
                                                         gene_to_uniprot_dict, 
                                                         cmap_path,
                                                         miss_prob_dict,
                                                         fragment=fragment,
                                                         alpha=alpha_pos,
                                                         num_iteration=num_iteration)

                # Go back to protein pos
                f_pos_result["Pos"] = f_pos_result["Pos"] + 1400 * (fragment-1)
                
                # Update community number according considering previous fragments
                f_community = len(f_pos_result["Community"][pd.notnull(f_pos_result["Community"])].unique())
                f_pos_result["Community"] = f_pos_result["Community"] + prot_community
                prot_community += f_community

                # Save fragments result
                f_pos_result_lst.append(f_pos_result)
                f_result_gene_lst.append(f_result_gene)


            result_pos = pd.concat(f_pos_result_lst).sort_values("pval", ascending=False).reset_index(drop=True)

            if pos_result is not None:
                result_pos_lst.append(pos_result)
            result_gene_lst.append(result_gene)

        ########### XXXXX

        











        # # Get contact map
        # try:
        #     # Here change to different F for bigger proteins
        #     cmap_gene = np.load(f"{cmap_path}/{gene_to_uniprot_dict[gene]}-F1.npy")                             #### Change the F if we want to include fragmented proteins

        #      # Run clustering
        #     try:
        #         pos_result, result_gene = pseudo_hotmaps(gene, 
        #                                         data,
        #                                         gene_to_uniprot_dict, 
        #                                         cmap_gene,
        #                                         miss_prob_dict,
        #                                         alpha=alpha_pos,
        #                                         num_iteration=num_iteration)
        #         result_gene_lst.append(result_gene)
                
        #         if pos_result is not None:
        #             result_pos_lst.append(pos_result)
         
        #     # Can't run clustering                                                                >>>> Should raise a better exception to capture a more specific error
        #     except:
        #         result_gene = pd.DataFrame({"Gene" : gene,
        #                                 "No_frag" : np.nan,
        #                                 "Mut_count" : np.nan,
        #                                 "Max_mut_pos" : np.nan,
        #                                 "Structure_max_pos" : np.nan,
        #                                 "Status" : "Not_processed"},
        #                                 index=[0])
        #         result_gene_lst.append(result_gene)
            
        # # Cant load cmap
        # except:
        #     result_gene = pd.DataFrame({"Gene" : gene,
        #                            "No_frag" : np.nan,
        #                            "Mut_count" : np.nan,
        #                            "Max_mut_pos" : np.nan,
        #                            "Structure_max_pos" : np.nan,
        #                            "Status" : "Cmap_not_found"},
        #                             index=[0])
        #     result_gene_lst.append(result_gene)

    
    ## Save 
    print("\nSaving..")
    result_gene = pd.concat(result_gene_lst)

    # Plot 0                                                                                               ### TO REMOVE
    status_count = result_gene.groupby("Status").count().reset_index()[["Status", "Gene"]]
    status_count = status_count.rename(columns={"Gene" : "Count"})

    plt.figure(figsize=(12, 4))
    ax = sns.barplot(x='Status', y='Count', data=status_count)
    ax.bar_label(ax.containers[0])
    plt.savefig(f"{output_dir}/{cancer_type}.{cohort}.processing_summary.png", bbox_inches="tight", dpi=300)
    plt.clf()

    if len(result_pos_lst) == 0:
        print(f"Did not processed any genes\n")
        result_gene.to_csv(f"{output_dir}/{cancer_type}.{cohort}.3d_clustering_genes.csv", index=False)
    else:
        result_pos = pd.concat(result_pos_lst)
        result_pos.to_csv(f"{output_dir}/{cancer_type}.{cohort}.3d_clustering_pos.csv", index=False)

        # Get gene global pval, qval, and clustering annotations
        result_gene = get_final_gene_result(result_pos, result_gene, alpha_gene)
        result_gene.to_csv(f"{output_dir}/{cancer_type}.{cohort}.3d_clustering_genes.csv", index=False)


if __name__ == "__main__":
    main()
 