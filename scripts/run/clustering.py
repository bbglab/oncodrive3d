"""
Contains functions to perform the 3D clustering of missense mutations.
"""

import multiprocessing
import os

import daiquiri
import networkx.algorithms.community as nx_comm
import numpy as np
import pandas as pd

from scripts import __logger_name__
from scripts.run.communities import get_community_index_nx, get_network
from scripts.run.miss_mut_prob import get_unif_gene_miss_prob
from scripts.run.score_and_simulations import (get_anomaly_score,
                                               get_sim_anomaly_score,
                                               recompute_inf_score)
from scripts.run.utils import add_info

logger = daiquiri.getLogger(__logger_name__ + ".run.clustering")


def process_mapping_issue(issue_ix, 
                          mut_gene_df, 
                          result_gene_df, 
                          gene, uniprot_id, 
                          af_f, 
                          transcript_status, 
                          thr, 
                          issue_type):
    """
    Check if there are mutations not in the structure (mut pos exceed lenght of 
    the structure protein sequence) or mutations with mismatches between WT AA 
    between mut and structure protein sequence. If the ratio of mutations do 
    not exceed threshold, filter out the specific mutations, else filter out 
    all mutations of that gene.
    """

    if issue_type == "Mut_not_in_structure":
        logger_txt="mut not in the structure"
        df_col = "Ratio_not_in_structure"
    elif issue_type == "WT_mismatch":
        logger_txt="mut with ref-structure WT AA mismatch"
        df_col = "Ratio_WT_mismatch"
    else:
        logger.warning(f"'{issue_type}' is not a valid issue type, please select 'Mut_not_in_structure' or 'Ratio_not_in_structure': Skipping processing mapping issue..")
        filter_gene = False
        
        return filter_gene, result_gene_df, mut_gene_df
        
    ratio_issue = sum(issue_ix) / len(mut_gene_df)
    logger_out = f"Detected {sum(issue_ix)} ({ratio_issue*100:.1f}%) {logger_txt} of {gene} ({uniprot_id}-F{af_f}, transcript status = {transcript_status}): "
    result_gene_df[df_col] = ratio_issue
    
    if ratio_issue > thr:
        result_gene_df["Status"] = issue_type
        if transcript_status == "Match":
            logger.warning(logger_out + "Filtering the gene")
        else:
            logger.debug(logger_out + "Filtering the gene..")
        filter_gene = True
        
        return filter_gene, result_gene_df, None
    
    else:
        if transcript_status == "Match":
            logger.warning(logger_out + "Filtering the mutations..")
        else:
            logger.debug(logger_out + "Filtering the mutations..")    
        mut_gene_df = mut_gene_df[~issue_ix]
        filter_gene = False
        
        return filter_gene, result_gene_df, mut_gene_df
    
    
def clustering_3d(gene, 
                  uniprot_id,
                  mut_gene_df,                                      
                  cmap_path,
                  miss_prob_dict,
                  seq_gene,
                  af_f,
                  alpha=0.01,
                  num_iteration=10000,
                  cmap_prob_thr=0.5,
                  seed=None,
                  pae_path=None,
                  thr_mapping_issue=0.1):
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
    transcript_id_input = mut_gene_df.Transcript_ID.iloc[0]
    transcript_id_o3d = mut_gene_df.O3D_transcript_ID.iloc[0]
    transcript_status = mut_gene_df.Transcript_status.iloc[0]
    result_gene_df = pd.DataFrame({"Gene" : gene,
                                   "Uniprot_ID" : uniprot_id,
                                   "F" : af_f,                           
                                   "Mut_in_gene" : mut_count,
                                   "Ratio_not_in_structure" : 0,
                                   "Ratio_WT_mismatch" : 0,
                                   "Mut_zero_mut_prob" : 0,
                                   "Pos_zero_mut_prob" : np.nan,
                                   "Transcript_ID" : transcript_id_input,
                                   "O3D_transcript_ID" : transcript_id_o3d,
                                   "Transcript_status" : transcript_status,
                                   "Status" : np.nan}, 
                                    index=[1])


    # Check if there is a mutation that is not in the structure      
    if max(mut_gene_df.Pos) > len(seq_gene):
        not_in_structure_ix = mut_gene_df.Pos > len(seq_gene)
        filter_gene, result_gene_df, mut_gene_df = process_mapping_issue(not_in_structure_ix, 
                                                                         mut_gene_df, 
                                                                         result_gene_df, 
                                                                         gene, 
                                                                         uniprot_id, 
                                                                         af_f, 
                                                                         transcript_status, 
                                                                         thr_mapping_issue,
                                                                         issue_type="Mut_not_in_structure")
        if filter_gene:
            return None, result_gene_df
            
    # Check for mismatch between WT reference and WT structure 
    wt_mismatch_ix = mut_gene_df.apply(lambda x: seq_gene[x.Pos-1] != x.WT, axis=1)
    if sum(wt_mismatch_ix) > 0:
        filter_gene, result_gene_df, mut_gene_df = process_mapping_issue(wt_mismatch_ix, 
                                                                         mut_gene_df, 
                                                                         result_gene_df, 
                                                                         gene, 
                                                                         uniprot_id, 
                                                                         af_f, 
                                                                         transcript_status, 
                                                                         thr_mapping_issue,
                                                                         issue_type="WT_mismatch")
        if filter_gene:
            return None, result_gene_df

    # Load cmap
    cmap_complete_path = f"{cmap_path}/{uniprot_id}-F{af_f}.npy"
    if os.path.isfile(cmap_complete_path):
        cmap = np.load(cmap_complete_path) 
        cmap = cmap > cmap_prob_thr
        cmap = cmap.astype(int)
    else:
        result_gene_df["Status"] = "Cmap_not_found"
        return None, result_gene_df
    
    # Load PAE
    pae_complete_path = f"{pae_path}/{uniprot_id}-F{af_f}-predicted_aligned_error.npy"
    if os.path.isfile(pae_complete_path):
        pae = np.load(pae_complete_path) 
    else:
        pae = None
    

    ## Get expected local myssense mutation density

    # Probability that each residue can be hit by a missense mut
    if miss_prob_dict is not None:
        gene_miss_prob = np.array(miss_prob_dict[f"{uniprot_id}-F{af_f}"])
    else:
        gene_miss_prob = get_unif_gene_miss_prob(size=len(cmap))

    # Filter out genes whose missense prob vec include any NA
    if np.any(np.isnan(gene_miss_prob)):
        result_gene_df["Status"] = "NA_miss_prob"
        return None, result_gene_df
    
    # Filter out genes with a mutation in a residue having zero prob to mutate
    pos_vec = np.unique(mut_gene_df["Pos"].values)
    pos_prob_vec = np.array(gene_miss_prob)[pos_vec-1]
    if (pos_prob_vec == 0).any():
        result_gene_df["Status"] = "Mut_with_zero_prob"
        pos_prob_vec = np.array(gene_miss_prob)[pos_vec-1]
        pos_zero_prob = list(pos_vec[pos_prob_vec == 0])
        mut_zero_prob_ix = mut_gene_df["Pos"].isin(pos_zero_prob)
        mut_zero_prob_count = sum(mut_zero_prob_ix)
        result_gene_df["Mut_zero_mut_prob"] = mut_zero_prob_count
        result_gene_df["Pos_zero_mut_prob"] = str(pos_zero_prob)
        ratio_zero_prob = mut_zero_prob_count / len(mut_gene_df)
        logger_out = f"Detected {mut_zero_prob_count} ({ratio_zero_prob*100:.1f}%) mut in {len(pos_zero_prob)} pos {pos_zero_prob} with zero mut prob in {gene} ({uniprot_id}-F{af_f}, transcript status = {transcript_status}): "
        result_gene_df["Ratio_mut_zero_prob"] = ratio_zero_prob 
        
        if ratio_zero_prob > thr_mapping_issue:
            result_gene_df["Status"] = "Mut_with_zero_prob"
            if transcript_status == "Match":
                logger.warning(logger_out + "Filtering the gene..")
            else:
                logger.debug(logger_out + "Filtering the gene..")
            return None, result_gene_df
        else:
            if transcript_status == "Match":
                logger.warning(logger_out + "Filtering the mutations..")
            else:
                logger.debug(logger_out + "Filtering the mutations..")    
            mut_gene_df = mut_gene_df[~mut_zero_prob_ix]

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
    sim_anomaly = get_sim_anomaly_score(len(mut_gene_df), 
                                        cmap, 
                                        gene_miss_prob,
                                        vol_missense_mut_prob,                                                                               
                                        num_iteration=num_iteration,
                                        seed=seed) 

    # Get ranked observed score (loglik+_LFC) 
    no_mut_pos = len(result_pos_df)
    sim_anomaly = sim_anomaly.iloc[:no_mut_pos,:].reset_index()
    result_pos_df["Score"] = get_anomaly_score(result_pos_df["Mut_in_vol"], 
                                                     len(mut_gene_df), 
                                                     vol_missense_mut_prob[result_pos_df["Pos"]-1])
    if np.isinf(result_pos_df.Score).any():
        logger.debug(f"Detected inf observed score in gene {gene} ({uniprot_id}-F{af_f}): Recomputing with higher precision..")
        result_pos_df = recompute_inf_score(result_pos_df, len(mut_gene_df), vol_missense_mut_prob[result_pos_df["Pos"]-1])
    
    mut_in_res = count.rename("Mut_in_res").reset_index().rename(columns={"index" : "Pos"})
    result_pos_df = mut_in_res.merge(result_pos_df, on = "Pos", how = "outer")                          
    result_pos_df = result_pos_df.sort_values("Score", ascending=False).reset_index(drop=True)


    ## Compute p-val and assign hits

    # Add to the simulated score of each iteration its standard deviation  
    # (makes the method more conservative, eg., avoid borderline cases)
    sim_anomaly.iloc[:,1:] = sim_anomaly.apply(lambda x: x[1:] + x[1:].std(), axis=1)

    # Ratio observed and simulated anomaly scores 
    # (used to break the tie in p-values gene sorting)
    result_pos_df["Score_obs_sim"] = sim_anomaly.apply(lambda x: result_pos_df["Score"].values[int(x["index"])] / np.mean(x[1:]), axis=1) 

    # Empirical p-val
    result_pos_df["pval"] = sim_anomaly.apply(lambda x: sum(x[1:] >= result_pos_df["Score"].values[int(x["index"])]) / len(x[1:]), axis=1)

    # Assign hits
    result_pos_df["C"] = [int(i) for i in result_pos_df["pval"] < alpha]              
    
    # Select extended significant hits
    pos_hits = result_pos_df[result_pos_df["C"] == 1].Pos
    neigh_pos_hits = list(set([pos for p in pos_hits.values for pos in list(np.where(cmap[p - 1])[0] + 1)]))
    pos_hits_extended = [pos for pos in result_pos_df.Pos if pos in neigh_pos_hits]
    result_pos_df["C_ext"] = result_pos_df.apply(lambda x: 1 if (x["C"] == 0) & (x["Pos"] in pos_hits_extended)
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
        clustered_mut = sum([pos in np.unique(np.concatenate([np.where(cmap[pos-1])[0]+1 for pos in pos_hits.values])) 
                             for pos in mut_gene_df.Pos])
    else:
        clustered_mut = 0
    result_pos_df["Rank"] = result_pos_df.index
    result_pos_df.insert(0, "Gene", gene)
    result_pos_df.insert(1, "Uniprot_ID", uniprot_id)
    result_pos_df.insert(2, "F", af_f)
    result_pos_df.insert(4, "Mut_in_gene", len(mut_gene_df))    
    result_pos_df = add_info(mut_gene_df, result_pos_df, cmap, pae)
    result_gene_df["Clust_res"] = len(pos_hits)
    result_gene_df["Clust_mut"] = clustered_mut
    result_gene_df["Status"] = "Processed"

    return result_pos_df, result_gene_df


def clustering_3d_mp(genes,
                     data,
                     cmap_path,
                     miss_prob_dict,
                     seq_df,
                     plddt_df,
                     num_process,
                     alpha=0.01,
                     num_iteration=10000,
                     cmap_prob_thr=0.5,
                     seed=None,
                     pae_path=None,
                     thr_mapping_issue=0.1):
    """
    Run the 3D-clustering algorithm in parallel on multiple genes.
    """
    
    result_gene_lst = []
    result_pos_lst = []
    
    for n, gene in enumerate(genes):
    
        mut_gene_df = data[data["Gene"] == gene]
        seq_df_gene = seq_df[seq_df["Gene"] == gene]
        uniprot_id = seq_df_gene['Uniprot_ID'].values[0]
        seq = seq_df_gene['Seq'].values[0]
        af_f = seq_df_gene['F'].values[0]
        
        # Add confidence to mut_gene_df
        plddt_df_gene_df = plddt_df[plddt_df["Uniprot_ID"] == uniprot_id].drop(columns=["Uniprot_ID"])
        mut_gene_df = mut_gene_df.merge(plddt_df_gene_df, on = ["Pos"], how = "left")

        pos_result, result_gene = clustering_3d(gene,
                                                uniprot_id, 
                                                mut_gene_df, 
                                                cmap_path,
                                                miss_prob_dict,
                                                seq_gene=seq,
                                                af_f=af_f,
                                                alpha=alpha,
                                                num_iteration=num_iteration,
                                                cmap_prob_thr=cmap_prob_thr,
                                                seed=seed,
                                                pae_path=pae_path,
                                                thr_mapping_issue=thr_mapping_issue)
        result_gene_lst.append(result_gene)
        if pos_result is not None:
            result_pos_lst.append(pos_result)
            
        # Monitor processing
        if n == 0:
            logger.debug(f"Process [{num_process+1}] starting..")
        elif n % 10 == 0:
            logger.debug(f"Process [{num_process+1}] completed [{n+1}/{len(genes)}] structures..")
        elif n+1 == len(genes):
            logger.debug(f"Process [{num_process+1}] completed!")

    return result_gene_lst, result_pos_lst


def clustering_3d_mp_wrapper(genes,
                             data,
                             cmap_path,
                             miss_prob_dict,
                             seq_df,
                             plddt_df,
                             num_cores,
                             alpha=0.01,
                             num_iteration=10000,
                             cmap_prob_thr=0.5,
                             seed=None,
                             pae_path=None,
                             thr_mapping_issue=0.1):
    """
    Wrapper function to run the 3D-clustering algorithm in parallel on multiple genes.
    """

    # Split the genes into chunks for each process
    chunk_size = int(len(genes) / num_cores) + 1
    chunks = [genes[i : i + chunk_size] for i in range(0, len(genes), chunk_size)]
    # num_cores = min(num_cores, len(chunks))
    
    # Create a pool of processes and run clustering in parallel
    with multiprocessing.Pool(processes = num_cores) as pool:
        
        logger.debug(f'Starting [{len(chunks)}] processes..')
        results = pool.starmap(clustering_3d_mp, [(chunk,
                                                   data[data["Gene"].isin(chunk)], 
                                                   cmap_path, 
                                                   miss_prob_dict, 
                                                   seq_df[seq_df["Gene"].isin(chunk)],
                                                   plddt_df[plddt_df["Uniprot_ID"].isin(seq_df.loc[seq_df["Gene"].isin(chunk), "Uniprot_ID"])],
                                                   n_process,
                                                   alpha, 
                                                   num_iteration, 
                                                   cmap_prob_thr, 
                                                   seed, 
                                                   pae_path,
                                                   thr_mapping_issue) 
                                                  for n_process, chunk in enumerate(chunks)])
        
    # Parse output
    result_pos_lst = [pd.concat(r[1]) for r in results if len(r[1]) > 0]
    if len(result_pos_lst) > 0: 
        result_pos = pd.concat(result_pos_lst)
    else:
        result_pos = None
    result_gene = pd.concat([pd.concat(r[0]) for r in results])
    
    return result_pos, result_gene