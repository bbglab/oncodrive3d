"""
Contains functions to perform the 3D clustering of missense mutations.
"""

import multiprocessing
import os
import json
import daiquiri
import networkx.algorithms.community as nx_comm
import numpy as np
import pandas as pd

from scripts import __logger_name__

from scripts.run.communities import get_community_index_nx, get_network
from scripts.run.pvalues import get_final_gene_result
from scripts.run.miss_mut_prob import get_miss_mut_prob_dict, mut_rate_vec_to_dict, get_unif_gene_miss_prob
from scripts.run.score_and_simulations import get_anomaly_score, get_sim_anomaly_score, recompute_inf_score
from scripts.run.utils import add_info, get_gene_entry, add_nan_clust_cols, parse_maf_input, sort_cols, empty_result_pos
from scripts.run.mutability import init_mutabilities_module

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
    result_gene_df[df_col] = np.round(ratio_issue, 3)
    
    
    # Do not filter neither gene neither mut
    if issue_type == "WT_mismatch" and thr == 1:
        logger.warning(logger_out + "Filtering of mismatching mutations disabled ('thr_mapping_issue' = 1)..")
        mut_gene_df.loc[issue_ix, "WT_mismatch"] = 1 
        filter_gene = False
        return filter_gene, result_gene_df, mut_gene_df

    else:
        # Filter gene
        if ratio_issue >= thr:
            result_gene_df["Status"] = issue_type
            
            if transcript_status == "Match":
                logger.warning(logger_out + "Filtering the gene")
            else:
                logger.debug(logger_out + "Filtering the gene..")
            filter_gene = True
            return filter_gene, result_gene_df, None
        
        # Filter mut
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
                  thr_mapping_issue=0.1,
                  sample_info=False):
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
    wt_mismatch_ix = mut_gene_df.apply(lambda x: bool(seq_gene[x.Pos-1] != x.WT), axis=1)
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
            clumps = get_community_index_nx(pos_hits, communities)

        else:
            # Assign cluster 0 to the only pos hit
            clumps = 0 
        meta_clusters = pd.DataFrame({"Pos" : pos_hits, "Clump" : clumps})
        result_pos_df = result_pos_df.merge(meta_clusters, how = "left", on = "Pos")
    else:
        result_pos_df["Clump"] = np.nan
    

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
    result_pos_df = add_info(mut_gene_df, result_pos_df, cmap, pae, sample_info)
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
                     thr_mapping_issue=0.1,
                     sample_info=False):
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
                                                thr_mapping_issue=thr_mapping_issue,
                                                sample_info=sample_info)
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
                             thr_mapping_issue=0.1,
                             sample_info=False):
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
                                                   thr_mapping_issue,
                                                   sample_info) 
                                                  for n_process, chunk in enumerate(chunks)])
        
    # Parse output
    result_pos_lst = [pd.concat(r[1]) for r in results if len(r[1]) > 0]
    if len(result_pos_lst) > 0: 
        result_pos = pd.concat(result_pos_lst)
    else:
        result_pos = None
    result_gene = pd.concat([pd.concat(r[0]) for r in results])
    
    return result_pos, result_gene


def run_clustering(input_path,
                    mut_profile_path,
                    mutability_config_path,
                    output_dir,
                    cmap_path,
                    seq_df_path,
                    plddt_path,
                    pae_path,
                    n_iterations,
                    alpha,
                    cmap_prob_thr,
                    cores,
                    seed,
                    verbose,
                    cancer_type,
                    cohort,
                    no_fragments,
                    only_processed,
                    thr_mapping_issue,
                    o3d_transcripts,
                    use_input_symbols,
                    mane,
                    sample_info):
    """
    Main function to lunch the 3D clustering analysis.
    """

    # Load
    # ====

    seq_df = pd.read_csv(seq_df_path, sep="\t")
    data, seq_df = parse_maf_input(input_path, 
                                    seq_df, 
                                    use_o3d_transcripts=o3d_transcripts,
                                    use_input_symbols=use_input_symbols, 
                                    mane=mane)

    if len(data) > 0:

        # Run
        # ===

        # Get genes with enough mut
        result_np_gene_lst = []
        genes = data.groupby("Gene").apply(len)
        genes_mut = genes[genes >= 2]
        genes_no_mut = genes[genes < 2].index

        if len(genes_no_mut) > 0:
            logger.debug(f"Detected [{len(genes_no_mut)}] genes without enough mutations: Skipping..")
            result_gene = pd.DataFrame({"Gene" : genes_no_mut,
                                        "Uniprot_ID" : np.nan,
                                        "F" : np.nan,
                                        "Mut_in_gene" : 1,
                                        "Ratio_not_in_structure" : np.nan,
                                        "Ratio_WT_mismatch" : np.nan,
                                        "Mut_zero_mut_prob" : np.nan,
                                        "Pos_zero_mut_prob" : np.nan,
                                        "Transcript_ID" : get_gene_entry(data, genes_no_mut, "Transcript_ID"),
                                        "O3D_transcript_ID" : get_gene_entry(data, genes_no_mut, "O3D_transcript_ID"),
                                        "Transcript_status" : get_gene_entry(data, genes_no_mut, "Transcript_status"),
                                        "Status" : "No_mut"})
            result_np_gene_lst.append(result_gene)

        # Seq df for metadata info
        metadata_cols = [col for col in ["Gene", "HGNC_ID", "Ens_Gene_ID", "Ens_Transcr_ID", "Refseq_prot", "Uniprot_ID", "F"] if col in seq_df.columns]
        metadata_mapping_cols = [col for col in ["Seq", "Chr", "Reverse_strand", "Exons_coord", "Seq_dna", "Tri_context", "Reference_info"] if col in seq_df.columns]
        seq_df_all = seq_df[seq_df["Gene"].isin(genes.index)].copy()

        # Get genes with corresponding Uniprot-ID mapping        
        gene_to_uniprot_dict = {gene : uni_id for gene, uni_id in seq_df[["Gene", "Uniprot_ID"]].drop_duplicates().values}
        genes_to_process = [gene for gene in genes_mut.index if gene in gene_to_uniprot_dict.keys()]
        seq_df = seq_df[seq_df["Gene"].isin(genes_to_process)].reset_index(drop=True)
        genes_no_mapping = genes[[gene in genes_mut.index and gene not in gene_to_uniprot_dict.keys() for gene in genes.index]]
        if len(genes_no_mapping) > 0:
            logger.debug(f"Detected [{len(genes_no_mapping)}] genes without IDs mapping: Skipping..")
            result_gene = pd.DataFrame({"Gene" : genes_no_mapping.index,
                                        "Uniprot_ID" : np.nan,
                                        "F" : np.nan,
                                        "Mut_in_gene" : genes_no_mapping.values,
                                        "Ratio_not_in_structure" : np.nan,
                                        "Ratio_WT_mismatch" : np.nan,
                                        "Mut_zero_mut_prob" : np.nan,
                                        "Pos_zero_mut_prob" : np.nan,
                                        "Transcript_ID" : get_gene_entry(data, genes_no_mapping.index, "Transcript_ID"),
                                        "O3D_transcript_ID" : get_gene_entry(data, genes_no_mapping.index, "O3D_transcript_ID"),
                                        "Transcript_status" : get_gene_entry(data, genes_no_mapping.index, "Transcript_status"),
                                        "Status" : "No_ID_mapping"})
            result_np_gene_lst.append(result_gene)
        
        # Filter on fragmented (AF-F) genes
        if no_fragments:
            # Return the fragmented genes as non processed output
            genes_frag = seq_df[seq_df.F.str.extract(r'(\d+)', expand=False).astype(int) > 1]
            genes_frag = genes_frag.Gene.reset_index(drop=True).values
            genes_frag_mut = genes_mut[[gene in genes_frag for gene in genes_mut.index]]
            genes_frag = genes_frag_mut.index.values
            if len(genes_frag) > 0:
                logger.debug(f"Detected [{len(genes_frag)}] fragmented genes with disabled fragments processing: Skipping..")
                result_gene = pd.DataFrame({"Gene" : genes_frag,
                                            "Uniprot_ID" : np.nan, 
                                            "F" : np.nan,
                                            "Mut_in_gene" : genes_frag_mut.values,
                                            "Ratio_not_in_structure" : np.nan,
                                            "Ratio_WT_mismatch" : np.nan,
                                            "Mut_zero_mut_prob" : np.nan,
                                            "Pos_zero_mut_prob" : np.nan,
                                            "Transcript_ID" : get_gene_entry(data, genes_frag, "Transcript_ID"),
                                            "O3D_transcript_ID" : get_gene_entry(data, genes_frag, "O3D_transcript_ID"),
                                            "Transcript_status" : get_gene_entry(data, genes_frag, "Transcript_status"),
                                            "Status" : "Fragmented"})
                result_np_gene_lst.append(result_gene)
                # Filter out from genes to process and seq df
                genes_to_process = [gene for gene in genes_to_process if gene not in genes_frag]
                seq_df = seq_df[seq_df["Gene"].isin(genes_to_process)].reset_index(drop=True)
                
        # Filter on start-loss mutations
        start_mut_ix = data["Pos"] == 1
        start_mut = sum(start_mut_ix)
        if start_mut > 0:
            genes_start_mut = list(data[start_mut_ix].Gene.unique())
            data = data[~start_mut_ix]
            logger.warning(f"Detected {start_mut} start-loss mutations in {len(genes_start_mut)} genes {genes_start_mut}: Filtering mutations..")


        # Missense mut prob
        # =================
        
        # Using mutabilities if provided
        if mutability_config_path is not None:
            logger.info("Computing missense mut probabilities using mutabilities..")
            mutab_config = json.load(open(mutability_config_path, encoding="utf-8"))
            logger.debug("Init mutabilities module..")
            init_mutabilities_module(mutab_config)
            seq_df = seq_df[seq_df["Reference_info"] == 1]   
            seq_df['Exons_coord'] = seq_df['Exons_coord'].apply(eval)  
            genes_to_process = [gene for gene in genes_to_process if gene in seq_df["Gene"].unique()]
            genes_not_mutability = [gene for gene in genes_to_process if gene not in seq_df["Gene"].unique()]
            logger.debug("Computing probabilities..")
            miss_prob_dict = get_miss_mut_prob_dict(mut_rate_dict=None, seq_df=seq_df,
                                                    mutability=True, mutability_config=mutab_config)
            
            if len(genes_not_mutability) > 0:   
                logger.debug(f"Detected [{len(genes_not_mutability)}] genes without mutability information: Skipping..")
                result_gene = pd.DataFrame({"Gene" : genes_not_mutability,
                                            "Uniprot_ID" : np.nan,
                                            "F" : np.nan,
                                            "Mut_in_gene" : np.nan,
                                            "Ratio_not_in_structure" : np.nan,
                                            "Ratio_WT_mismatch" : np.nan,
                                            "Mut_zero_mut_prob" : np.nan,
                                            "Pos_zero_mut_prob" : np.nan,
                                            "Transcript_ID" : get_gene_entry(data, genes_not_mutability, "Transcript_ID"),
                                            "O3D_transcript_ID" : get_gene_entry(data, genes_not_mutability, "O3D_transcript_ID"),
                                            "Transcript_status" : get_gene_entry(data, genes_not_mutability, "Transcript_status"),
                                            "Status" : "No_mutability"})
                result_np_gene_lst.append(result_gene)
                
        # Using mutational profiles
        elif mut_profile_path is not None:
            # Compute dict from mut profile of the cohort and dna sequences
            mut_profile = json.load(open(mut_profile_path, encoding="utf-8"))
            logger.info("Computing missense mut probabilities..")
            if not isinstance(mut_profile, dict):
                mut_profile = mut_rate_vec_to_dict(mut_profile)
            miss_prob_dict = get_miss_mut_prob_dict(mut_rate_dict=mut_profile, seq_df=seq_df)
        else:
            logger.warning("Mutation profile not provided: Uniform distribution will be used for scoring and simulations.")
            miss_prob_dict = None


        # Run 3D-clustering
        # =================
        
        if len(result_np_gene_lst):
            result_np_gene = pd.concat(result_np_gene_lst)
            result_np_gene["Uniprot_ID"] = [gene_to_uniprot_dict[gene] if gene in gene_to_uniprot_dict.keys() else np.nan for gene in result_np_gene["Gene"].values]
        if len(genes_to_process) > 0:
            logger.info(f"Performing 3D-clustering on [{len(seq_df)}] proteins..")
            seq_df = seq_df[["Gene", "Uniprot_ID", "F", "Seq"]]
            plddt_df = pd.read_csv(plddt_path, sep="\t", usecols=["Pos", "Confidence", "Uniprot_ID"], dtype={"Pos" : np.int32,
                                                                                                                "Confidence" : np.float32, 
                                                                                                                "Uniprot_ID" : "object"})  
            
            result_pos, result_gene = clustering_3d_mp_wrapper(genes=genes_to_process,
                                                                data=data,
                                                                cmap_path=cmap_path,
                                                                miss_prob_dict=miss_prob_dict,
                                                                seq_df=seq_df,
                                                                plddt_df=plddt_df,
                                                                num_cores=cores,
                                                                alpha=alpha,
                                                                num_iteration=n_iterations,
                                                                cmap_prob_thr=cmap_prob_thr,
                                                                seed=seed,
                                                                pae_path=pae_path,
                                                                thr_mapping_issue=thr_mapping_issue,
                                                                sample_info=sample_info)
            if result_np_gene_lst:
                result_gene = pd.concat((result_gene, result_np_gene))
        else:
            result_gene = result_np_gene
            result_pos = None


        # Save
        #=====
        
        os.makedirs(output_dir, exist_ok=True)
        result_gene["Cancer"] = cancer_type
        result_gene["Cohort"] = cohort
        output_path_pos = os.path.join(output_dir, f"{cohort}.3d_clustering_pos.csv")
        output_path_genes = os.path.join(output_dir, f"{cohort}.3d_clustering_genes.csv")
        
        # Save processed seq_df and input files
        seq_df_output = os.path.join(output_dir, f"{cohort}.seq_df.processed.tsv")
        input_mut_output = os.path.join(output_dir, f"{cohort}.mutations.processed.tsv")
        input_prob_output = os.path.join(output_dir, f"{cohort}.miss_prob.processed.json")
        logger.info(f"Saving {seq_df_output}")
        seq_df_all[metadata_cols + metadata_mapping_cols].to_csv(seq_df_output, sep="\t", index=False)
        logger.info(f"Saving {input_mut_output}")
        data.to_csv(input_mut_output, sep="\t", index=False)
        logger.info(f"Saving {input_prob_output}")
        with open(input_prob_output, "w") as json_file:
            json.dump(miss_prob_dict, json_file)
        
        # Add extra metadata
        result_gene = result_gene.merge(seq_df_all[metadata_cols], on=["Gene", "Uniprot_ID"], how="left")

        if only_processed:
            result_gene = result_gene[result_gene["Status"] == "Processed"]

        if result_pos is None:
            # Save gene-level result and empty res-level result
            logger.warning("Did not processed any genes!")
            result_gene = add_nan_clust_cols(result_gene)
            result_gene = sort_cols(result_gene)
            if not sample_info:
                result_gene.drop(columns=[col for col in ['Tot_samples', 
                                                            'Samples_in_top_vol', 
                                                            'Samples_in_top_cl_vol'] if col in result_gene.columns], inplace=True)   
            if no_fragments:
                result_gene = result_gene.drop(columns=[col for col in ["F", "Mut_in_top_F", "Top_F"] if col in result_gene.columns])
            empty_result_pos(sample_info).to_csv(output_path_pos, index=False)
            result_gene.to_csv(output_path_genes, index=False)

            logger.info(f"Saving (empty) {output_path_pos}")
            logger.info(f"Saving {output_path_genes}")

        else:
            # Save res-level result
            result_pos["Cancer"] = cancer_type
            result_pos["Cohort"] = cohort
            if not sample_info:
                result_pos.drop(columns=[col for col in ['Tot_samples', 
                                                            'Samples_in_vol', 
                                                            'Samples_in_cl_vol'] if col in result_gene.columns], inplace=True) 
            result_pos = result_pos.sort_values(["Gene", "pval", "Score_obs_sim"], ascending=[True, True, False]).reset_index(drop=True)
            result_pos.to_csv(output_path_pos, index=False)

            # Get gene global pval, qval, and clustering annotations and save gene-level result
            result_gene = get_final_gene_result(result_pos, result_gene, alpha, sample_info)
            result_gene = sort_cols(result_gene) 
            if not sample_info:
                result_gene.drop(columns=[col for col in ['Tot_samples', 
                                                          'Samples_in_top_vol', 
                                                          'Samples_in_top_cl_vol'] if col in result_gene.columns], inplace=True)   
            if no_fragments:
                result_gene.drop(columns=[col for col in ["F", "Mut_in_top_F", "Top_F"] if col in result_gene.columns], inplace=True)
            with np.printoptions(linewidth=10000):
                result_gene.to_csv(output_path_genes, index=False)

            logger.info(f"Saving {output_path_pos}")
            logger.info(f"Saving {output_path_genes}")

        logger.info("3D-clustering analysis completed!")

    else:
        logger.warning("No missense mutations were found in the input MAF. Consider checking your data: the field 'Variant_Classification' should include either 'Missense_Mutation' or 'missense_variant'")