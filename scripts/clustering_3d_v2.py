""" 
The module includes the functions and the script to run 
an HotMAPs-inspired method that perform 3D-clustering 
of mutations using simulations, rank based comparison and
the canonical predicted structure stored in AlphaFold db.

###################################### EXAMPLE USAGE ############################################

time python3 clustering_3d.py -i ../tests/input/HARTWIG_WGS_PANCREAS_2020.in.maf -o ../tests/output/ \
-P ../tests/input/HARTWIG_WGS_PANCREAS_2020.mutrate.json -H 0


python3 clustering_3d.py -i ../../../evaluation/datasets/input/maf/PCAWG_WGS_COLORECT_ADENOCA.in.maf \         
-o ../../../evaluation/output \
-p ../../../evaluation/datasets/input/mut_profile/prior.weighted.normalized.PCAWG_WGS_COLORECT_ADENOCA.json \
-s ../datasets/seq_for_mut_prob.csv \
-c ../datasets/cmaps/ \
-u ../datasets/af_uniprot_to_gene_id.json \
-n 10000 \
-H 0 \
-t COREAD \
-C PCAWG_WGS_COLORECT_ADENOCA


#################################################################################################
"""


import os
import json
import numpy as np
# import seaborn as sns
# import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp
from progressbar import progressbar
import argparse
import networkx as nx
import networkx.algorithms.community as nx_comm
from scipy.stats import rankdata
from utils.utils import parse_maf_input, uniprot_to_hugo, get_pos_fragments
from utils.miss_mut_prob import mut_rate_vec_to_dict, get_miss_mut_prob_dict



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


## Process p-val and gene results

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
    result_gene["Max_mut_in_vol"] = -result_gene["Max_mut_in_vol"]
    result_gene = result_gene.sort_values(["pval", "Max_mut_in_vol", "Mut_in_gene"], ascending=True).reset_index(drop=True)
    result_gene["Max_mut_in_vol"] = -result_gene["Max_mut_in_vol"]
    
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
    result_gene.insert(4, "C_gene", gene_label)
    result_gene.insert(8, "Max_mut_in_vol", result_gene.pop("Mut_in_vol"))
    result_gene.insert(1, "Uniprot_ID", result_gene.pop("Uniprot_ID"))
    result_gene = sort_by_pval_and_mut(result_gene)
    
    return result_gene


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
                   uniprot_id,
                   mut_gene_df,                                      
                   cmap_path,
                   miss_prob_dict,
                   fragment=1,
                   alpha=0.01,
                   num_iteration=10000,
                   hits_only=True):
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
                                   "Uniprot_ID" : uniprot_id,
                                   "No_frag" : fragment,                           
                                   "Mut_in_gene" : mut_count,
                                   "Max_mut_pos" : np.nan,
                                   "Structure_max_pos" : np.nan,
                                   "Status" : np.nan}, 
                                    index=[1])
    
    # Load cmap
    cmap_complete_path = f"{cmap_path}/{uniprot_id}-F{fragment}.npy"
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
    gene_miss_prob = np.array(miss_prob_dict[f"{uniprot_id}-F{fragment}"])

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
        

    ## Output
    result_pos_df.insert(0, "Gene", gene)
    result_pos_df.insert(1, "Uniprot_ID", uniprot_id)
    result_pos_df.insert(2, "F", fragment)
    result_pos_df.insert(4, "Mut_in_gene", mut_count)
    result_gene_df["Status"] = "Processed"

    # Keep only positions in clusters           
    if hits_only:
        result_pos_df = result_pos_df[result_pos_df["C"] == 1]

    return result_pos_df, result_gene_df


def main():
    """
    Wrapper function.
    """

    ## Initialize
    version = "2023_v1.3"    # LAST CHANGE: obtain HUGO mapping from seq_df

    # Parser
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)

    parser.add_argument("-i", "--input_maf", help="Path of the maf file used as input", type=str, required=True)
    parser.add_argument("-o", "--output_dir", help="Path to output directory", type=str, required=True)

    group.add_argument("-p", "--mut_profile", help="Path to the mut profile (list of 96 floats) of the cohort (json)", type=str)
    group.add_argument("-P", "--miss_mut_prob", help="Path to the dict of missense mut prob of each protein based on mut profile of the cohort (json)",  
                       type=str)
    parser.add_argument("-s", "--seq_df", 
                        help="Path to the dataframe including DNA and protein seq of all gene/proteins (all AF predicted ones)", 
                        type=str, 
                        default="../datasets/seq_for_mut_prob.csv")       

    parser.add_argument("-c", "--cmap_path", help="Path to the directory containting the contact map of each protein", type=str, 
                        default="../datasets/cmaps/")

    parser.add_argument("-n", "--n_iterations", help="Number of densities to be simulated", type=int, default=10000)
    parser.add_argument("-a", "--alpha_level_res", help="Significant threshold for the p-value of protein residues", type=float, default=0.01)
    parser.add_argument("-A", "--alpha_level_gene", help="Significant threshold for the global q-value of the genes", type=float, default=0.05)
    parser.add_argument("-H", "--hits_only", help="if 1 returns only positions in clusters, if 0 returns all", type=int, default=1)
    
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
    hits_only = args.hits_only

    print(f"Starting 3D-clustering [{version}]..\n")
    print(f"Iterations: {num_iteration}")
    print(f"Significant level (position): {alpha_pos}")
    print(f"Global significant level (gene): {alpha_gene}")
    print(f"Cohorts: {cohort}")
    print(f"Cancer type: {cancer_type}")
    print(f"Output directory: {output_dir}")


    ## Load input and df of DNA sequences

    # MAF input
    data = parse_maf_input(maf_input_path)

    # Seq df for missense mut prob
    seq_df = pd.read_csv(seq_df_path)


    ## Run

    result_pos_lst = []
    result_gene_lst = []

    # Get genes with enough mut
    genes = data.groupby("Gene").apply(lambda x: len(x))
    genes_mut = genes[genes >= 2].index
    genes_no_mut = genes[genes < 2].index
    result_gene = pd.DataFrame({"Gene" : genes_no_mut,
                                "Uniprot_ID" : np.nan,
                                "No_frag" : np.nan,
                                "Mut_in_gene" : 1,
                                "Max_mut_pos" : np.nan,
                                "Structure_max_pos" : np.nan,
                                "Status" : "No_mut"})
    result_gene_lst.append(result_gene)

    # Get genes with corresponding Uniprot-ID mapping
    gene_to_uniprot_dict = {gene : uni_id for gene, uni_id in seq_df[["Gene", "Uniprot_ID"]].drop_duplicates().values}
    genes_mapped = [gene for gene in genes_mut if gene in gene_to_uniprot_dict.keys()]
    seq_df = seq_df[[gene in genes_mapped for gene in seq_df["Gene"]]].reset_index(drop=True)
    genes_no_mapping = genes[[gene in genes_mut and gene not in gene_to_uniprot_dict.keys() for gene in genes.index]]
    result_gene = pd.DataFrame({"Gene" : genes_no_mapping.index,
                                "Uniprot_ID" : np.nan,
                                "No_frag" : np.nan,
                                "Mut_in_gene" : genes_no_mapping.values,
                                "Max_mut_pos" : np.nan,
                                "Structure_max_pos" : np.nan,
                                "Status" : "No_ID_mapping"})
    result_gene_lst.append(result_gene)

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
    print("Performing 3D clustering on genes with enough mutations..")
    for gene in progressbar(genes_mapped):
        mut_gene_df = data[data["Gene"] == gene]
        uniprot_id = gene_to_uniprot_dict[gene]   

        # If there is a single fragment
        if seq_df[seq_df["Gene"] == gene].F.max() == 1:      
            
            # Run clustering
            try:
                pos_result, result_gene = pseudo_hotmaps(gene,
                                                         uniprot_id, 
                                                         mut_gene_df, 
                                                         cmap_path,
                                                         miss_prob_dict,
                                                         fragment=1,
                                                         alpha=alpha_pos,
                                                         num_iteration=num_iteration,
                                                         hits_only=hits_only)
                result_gene_lst.append(result_gene)
                
                if pos_result is not None:
                    result_pos_lst.append(pos_result)
        
            # Can't run clustering                                                                >>>> Should raise a better exception to capture a more specific error
            except:
                result_gene = pd.DataFrame({"Gene" : gene,
                                            "Uniprot_ID" : uniprot_id,
                                            "No_frag" : np.nan,
                                            "Mut_in_gene" : np.nan,
                                            "Max_mut_pos" : np.nan,
                                            "Structure_max_pos" : np.nan,
                                            "Status" : "Not_processed"},
                                            index=[0])
                result_gene_lst.append(result_gene)

        # If the protein is fragmented
        else:

            mut_gene_df.insert(len(mut_gene_df.columns), "F", get_pos_fragments(mut_gene_df))     
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
                                                             uniprot_id,
                                                             mut_fragment_df,
                                                             cmap_path,
                                                             miss_prob_dict,
                                                             fragment=fragment,
                                                             alpha=alpha_pos,
                                                             num_iteration=num_iteration,
                                                             hits_only=hits_only)

                if f_pos_result is not None:

                    # Go back to protein pos
                    f_pos_result["Pos"] = f_pos_result["Pos"] + 1400 * (fragment-1)
                    
                    # Update community number according considering previous fragments
                    f_community = len(f_pos_result["Community"][pd.notnull(f_pos_result["Community"])].unique())
                    f_pos_result["Community"] = f_pos_result["Community"] + prot_community
                    prot_community += f_community

                    # Save fragments result
                    f_pos_result_lst.append(f_pos_result)
                f_result_gene_lst.append(f_result_gene)

            # Concat df for clustering of positions
            if len(f_pos_result_lst) > 0:
                pos_result = pd.concat(f_pos_result_lst)
            else:
                pos_result = None
            if pos_result is not None:
                result_pos_lst.append(pos_result)

            # Obtain a single df for the gene summary
            f_result_gene["No_frag"] = f_result_gene["No_frag"].max()
            f_result_gene["Mut_in_gene"] = f_result_gene["Mut_in_gene"].sum()
            result_gene = pd.DataFrame(f_result_gene.apply(lambda x: x.unique()[0] if len(x.dropna().unique()) < 2 else x.unique(), axis=0)).T
            result_gene_lst.append(result_gene)

    
    ## Save 
    print("\nSaving..")
    result_gene = pd.concat(result_gene_lst)
    result_gene["Cancer"] = cancer_type
    result_gene["Cohort"] = cohort

    # # Plot 0                                                                                               ### TO REMOVE
    # status_count = result_gene.groupby("Status").count().reset_index()[["Status", "Gene"]]
    # status_count = status_count.rename(columns={"Gene" : "Count"})

    # plt.figure(figsize=(12, 4))
    # ax = sns.barplot(x='Status', y='Count', data=status_count)
    # ax.bar_label(ax.containers[0])
    # plt.savefig(f"{output_dir}/{cancer_type}.{cohort}.processing_summary.png", bbox_inches="tight", dpi=300)
    # plt.clf()

    if len(result_pos_lst) == 0:
        print(f"Did not processed any genes\n")
        result_gene.to_csv(f"{output_dir}/{cancer_type}.{cohort}.3d_clustering_genes.csv", index=False)
    else:
        result_pos = pd.concat(result_pos_lst)
        result_pos["Cancer"] = cancer_type
        result_pos["Cohort"] = cohort
        result_pos.to_csv(f"{output_dir}/{cancer_type}.{cohort}.3d_clustering_pos.csv", index=False)

        # Get gene global pval, qval, and clustering annotations
        result_gene = get_final_gene_result(result_pos, result_gene, alpha_gene)
        result_gene.to_csv(f"{output_dir}/{cancer_type}.{cohort}.3d_clustering_genes.csv", index=False)


if __name__ == "__main__":
    main()
 