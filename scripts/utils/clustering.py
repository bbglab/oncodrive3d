"""
Contains functions to perform the 3D clustering of missense mutations.
"""

import os
import numpy as np
import pandas as pd
import networkx.algorithms.community as nx_comm
from utils.score_and_simulations import get_anomaly_score, get_sim_anomaly_score
from utils.communities import get_network, get_community_index_nx
from utils.utils import get_pos_fragments, get_samples_info, add_samples_info


def clustering_3d(gene, 
                   uniprot_id,
                   mut_gene_df,                                      
                   cmap_path,
                   miss_prob_dict,
                   alpha=0.01,
                   num_iteration=10000,
                   hits_only=True,
                   ext_hits=True):
    """
    Compute local density of missense mutations for a sphere of 10A around          
    each amino acid position of the selected gene product. It performed a 
    rank-based comparison between observed density and simulated ones in 
    absense of positive selection (cohort mut profile). Get an experimental 
    per-residue p-val for the local enrichment and a global p-val for the gene,
    which correspond to the minimum p-val across its positions.
    
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
    af_f = mut_gene_df.AF_F.unique()[0]
    result_gene_df = pd.DataFrame({"Gene" : gene,
                                   "Uniprot_ID" : uniprot_id,
                                   "F" : af_f,                           
                                   "Mut_in_gene" : mut_count,
                                   "Max_mut_pos" : np.nan,
                                   "Structure_max_pos" : np.nan,
                                   "Status" : np.nan}, 
                                    index=[1])

    # Load cmap
    cmap_complete_path = f"{cmap_path}/{uniprot_id}-F{af_f}.npy"
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


    # Samples info
    samples_info = get_samples_info(mut_gene_df, cmap)
    

    ## Get expected local myssense mutation density

    # Probability that each residue can be hit by a missense mut
    gene_miss_prob = np.array(miss_prob_dict[f"{uniprot_id}-F{af_f}"])

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
    result_pos_df["Obs_anomaly"] = get_anomaly_score(result_pos_df["Mut_in_vol"], mut_count, vol_missense_mut_prob[result_pos_df["Pos"]-1])
    mut_in_res = count.reset_index().rename(columns = {"Pos" : "Mut_in_res", "index" : "Pos"})       # (unnecessary but can be informative)
    result_pos_df = mut_in_res.merge(result_pos_df, on = "Pos", how = "outer")                       #     "        "        "        "     
    result_pos_df = result_pos_df.sort_values("Obs_anomaly", ascending=False).reset_index(drop=True)


    ## Compute p-val and assign hits

    # Add to the simulated score of each iteration its standard deviation  
    # (makes the method more conservative, eg., avoid borderline cases)
    sim_anomaly.iloc[:,1:] = sim_anomaly.apply(lambda x: x[1:] + x[1:].std(), axis=1)

    # Ratio observed and simulated anomaly scores 
    # (used to break the tie in p-values gene sorting)
    result_pos_df["Ratio_obs_sim"] = sim_anomaly.apply(lambda x: result_pos_df["Obs_anomaly"].values[int(x["index"])] / np.mean(x[1:]), axis=1) 
    #result_pos_df["Diff_obs_sim"] = sim_anomaly.apply(lambda x: (result_pos_df["Obs_anomaly"].values[int(x["index"])] - np.mean(x[1:])) / np.std(x[1:]), axis=1)  ##### TO REMOV ONE OF THEM

    # Empirical p-val
    result_pos_df["pval"] = sim_anomaly.apply(lambda x: sum(x[1:] >= result_pos_df["Obs_anomaly"].values[int(x["index"])]) / len(x[1:]), axis=1)

    # Assign hits
    result_pos_df["C"] = [int(i) for i in result_pos_df["pval"] < alpha]              
    
    # Select extended significant hits
    pos_hits = result_pos_df[result_pos_df["C"] == 1].Pos
    if ext_hits:
        anomaly_ratio_th = 0
    else:
        anomaly_ratio_th = 1
    neigh_pos_hits = list(set([pos for p in pos_hits.values for pos in list(np.where(cmap[p - 1])[0] + 1)]))
    pos_hits_extended = [pos for pos in result_pos_df.Pos if pos in neigh_pos_hits]
    result_pos_df["C_ext"] = result_pos_df.apply(lambda x: 1 if (x["C"] == 0) & (x["Pos"] in pos_hits_extended) & (x["Ratio_obs_sim"] > anomaly_ratio_th) 
                                                    else 0 if (x["C"] == 1) else np.nan, axis=1)
    result_pos_df["C"] = result_pos_df.apply(lambda x: 1 if (x["C"] == 1) | (x["C_ext"] == 1) else 0, axis=1)
    pos_hits = result_pos_df[result_pos_df["C"] == 1].Pos  

    
    ## Communities detection
    if len(pos_hits) > 0:
        if len(pos_hits) > 1:
            # Build network and perform detection
            G = get_network(pos_hits, mut_count_v, cmap)
            communities = nx_comm.label_propagation_communities(G)
            clusters = get_community_index_nx(pos_hits, communities)

        else:
            # Assign cluster 0 to the only pos hit
            clusters = 0 
        meta_clusters = pd.DataFrame({"Pos" : pos_hits, "Cluster" : clusters})
        result_pos_df = result_pos_df.merge(meta_clusters, how = "left", on = "Pos")
    else:
        result_pos_df["Cluster"] = np.nan
    

    ## Output
    if len(pos_hits) > 0:
        clustered_mut = sum([pos in np.unique(np.concatenate([np.where(cmap[pos-1])[0]+1 for pos in pos_hits.values])) for pos in mut_gene_df.Pos])
    else:
        clustered_mut = 0
    result_pos_df["Rank"] = result_pos_df.index
    result_pos_df.insert(0, "Gene", gene)
    result_pos_df.insert(1, "Uniprot_ID", uniprot_id)
    result_pos_df.insert(2, "F", af_f)
    result_pos_df.insert(4, "Mut_in_gene", mut_count)    
    result_pos_df = add_samples_info(mut_gene_df, result_pos_df, samples_info, cmap)
    result_gene_df["Clust_res"] = len(pos_hits)
    result_gene_df["Clust_mut"] = clustered_mut
    result_gene_df["Status"] = "Processed"

    # Keep only positions in clusters           
    if hits_only:
        result_pos_df = result_pos_df[result_pos_df["C"] == 1]

    return result_pos_df, result_gene_df