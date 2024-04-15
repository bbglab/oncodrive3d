"""
Module to generate a pandas dataframe including identifiers 
mapped to protein and DNA sequences.

The functions are used to extract all protein sequences of PDB 
structures in a given directory; use EMBOSS backtranseq back 
translate proteins sequences into DNA; generate a dataframe 
including HUGO symbol, Uniprot_ID, protein, and DNA sequences. 
This dataframe is required to get the probability of each residue 
to mutate (missense mutation) based on the mutation profile 
(mutation rate in 96 trinucleotide contexts) of the cohort. 
The per-residue missense mutation probability of each protein is 
then used to get the probability of a certain volume to be hit by 
a missense mutation.
"""


import json
import os
import re
import time

import ast
import daiquiri
import numpy as np
import multiprocessing
import pandas as pd
import requests
import sys
from tqdm import tqdm
from bgreference import hg38, mm39
from Bio.Seq import Seq
from urllib.request import urlretrieve


from scripts import __logger_name__
from scripts.datasets.utils import (download_single_file,
                                    get_af_id_from_pdb,
                                    get_pdb_path_list_from_dir,
                                    get_seq_from_pdb, 
                                    get_seq_similarity,
                                    translate_dna, 
                                    uniprot_to_hugo, 
                                    uniprot_to_hugo_pressed)

logger = daiquiri.getLogger(__logger_name__ + ".build.seq_for_mut_prob")


#===========
# Initialize
#===========

def initialize_seq_df(input_path, uniprot_to_gene_dict):
    """
    Parse any PDB structure from a given directory and create 
    a dataframe including HUGO symbol, Uniprot-ID, AF fragment,
    and protein sequence.
    """
    
    list_prot_path = get_pdb_path_list_from_dir(input_path)
    gene_lst = []
    uni_id_lst = []
    f_lst = []
    seq_lst = []

    for path_structure in tqdm(list_prot_path, total=len(list_prot_path), desc="Generating sequence df"):
        uniprot_id, f = get_af_id_from_pdb(path_structure).split("-F")
        if uniprot_id in uniprot_to_gene_dict:
            gene = uniprot_to_gene_dict[uniprot_id]
        else:
            gene = np.nan
        seq = "".join(list(get_seq_from_pdb(path_structure)))
        gene_lst.append(gene)
        uni_id_lst.append(uniprot_id)
        f_lst.append(f)
        seq_lst.append(seq)

    seq_df = pd.DataFrame({"Gene" : gene_lst, 
                           "Uniprot_ID" : uni_id_lst, 
                           "F" : f_lst, 
                           "Seq" : seq_lst
                           }).sort_values(["Gene", "F"])

    return seq_df


# def initialize_seq_df(path_to_pdb, uniprot_to_gene_dict):
#     """
#     Parse any PDB structure from a given directory and create 
#     a dataframe including HUGO symbol, Uniprot-ID, AF fragment,
#     and protein sequence.
#     """
    
#     # Get all PDB path in directory whose IDs can be mapped to HUGO symbols
#     list_prot_path = get_pdb_path_list_from_dir(path_to_pdb)
#     list_prot_path = [path for path in list_prot_path if get_af_id_from_pdb(path).split("-F")[0] in uniprot_to_gene_dict.keys()]
#     pdb_not_mapped = set([get_af_id_from_pdb(path).split("-F")[0] 
#                           for path in list_prot_path if get_af_id_from_pdb(path).split("-F")[0] not in uniprot_to_gene_dict.keys()])
#     if len(pdb_not_mapped) > 0:                                 
#         logger.warning(f"{len(pdb_not_mapped)} Uniprot-ID not found in the Uniprot-HUGO mapping dictionary.")

#     # Get Uniprot ID, HUGO, F and protein sequence of any PDB in dir
#     gene_lst = []
#     uni_id_lst = []
#     f_lst = []
#     seq_lst = []

#     for path_structure in tqdm(list_prot_path, total=len(list_prot_path), desc="Generating sequence df"):
#         uniprot_id, f = get_af_id_from_pdb(path_structure).split("-F")
#         gene = uniprot_to_gene_dict[uniprot_id]
#         seq = "".join(list(get_seq_from_pdb(path_structure)))
#         gene_lst.append(gene)
#         uni_id_lst.append(uniprot_id)
#         f_lst.append(f)
#         seq_lst.append(seq)

#     seq_df = pd.DataFrame({"Gene" : gene_lst, 
#                            "Uniprot_ID" : uni_id_lst, 
#                            "F" : f_lst, 
#                            "Seq" : seq_lst
#                            }).sort_values(["Gene", "F"])

#     return seq_df
    

# import shutil  ##  TO DELETE

# def initialize_seq_df_from_mane(path_to_datasets):
#     """
#     Parse any PDB structure from a given directory and create a 
#     dataframe including HUGO symbol, HGNC_ID, RefSeq_prot, MANE 
#     missing status, Uniprot-ID, AF fragment, and protein sequence.
#     Include only Uniprot IDs overlapping with MANE trascripts.
#     Mane missing status of 1 indicates that there might not be a perfect 
#     match between the Uniprot ID and the corresponding MANE transcript.
#     """

#     # Select structures
#     path_to_pdb = os.path.join(path_to_datasets, "pdb_structures")
#     uniprot_ids = os.listdir(path_to_pdb)
#     uniprot_ids = [uni_id.split("-")[1] for uni_id in list(set(uniprot_ids)) if ".pdb" in uni_id]

#     # Load MANE metadata
#     mane_summary_path = os.path.join(path_to_datasets, "mane_summary.txt.gz")
#     if not os.path.exists(mane_summary_path):
#         url_mane_summary = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.0/MANE.GRCh38.v1.0.summary.txt.gz"
#         _ = urlretrieve(url_mane_summary, f"{path_to_datasets}/mane_summary.txt.gz")
#     mane_to_af_path = os.path.join(path_to_datasets, "mane_refseq_prot_to_alphafold.csv")
#     mane_missing_path = os.path.join(path_to_datasets, "mane_missing.csv")
#     for path in [mane_to_af_path, mane_missing_path]:
#         if not os.path.exists(path):
#             logger.error(f"MANE metadata file {path} not found. Exiting..")
#     mane_to_af = pd.read_csv(mane_to_af_path)
#     mane_missing = pd.read_csv(mane_missing_path)
#     mane_summary = pd.read_csv(mane_summary_path, compression='gzip', sep="\t")

#     # Parse
#     mane_to_af = mane_to_af.rename(columns={"refseq_prot" : "RefSeq_prot",
#                                             "uniprot_accession" : "Uniprot_ID"})
#     mane_to_af["Mane_missing"] = 0
#     mane_to_af_missing = mane_missing.rename(columns={"refseq_prot" : "RefSeq_prot", 
#                                                       "uniprot_accession(s)" : "Uniprot_ID"})
#     mane_to_af_missing["Mane_missing"] = 1
#     mane_summary = mane_summary.rename(columns={"symbol" : "Gene"})
#     mane_to_af = pd.concat((mane_to_af[["RefSeq_prot", "Uniprot_ID", "Mane_missing"]], 
#                             mane_to_af_missing[["RefSeq_prot", "Uniprot_ID", "Mane_missing"]])).dropna(subset="Uniprot_ID")
#     id_mapping = mane_to_af[["RefSeq_prot", "Uniprot_ID", "Mane_missing"]].merge(
#         mane_summary[["RefSeq_prot", "Gene", "HGNC_ID"]], on="RefSeq_prot", how="inner")
#     id_mapping["Uniprot_ID"] = id_mapping.apply(lambda x: 
#                                                 select_uni_id(x.Uniprot_ID.split(";"), uniprot_ids) if len(x.Uniprot_ID.split(";")) > 1 
#                                                 else x.Uniprot_ID, axis=1)

#     # Select only available Uniprot IDs matching MANE
#     id_mapping = id_mapping.dropna(subset=["Uniprot_ID", "Gene"])
#     id_mapping = id_mapping[id_mapping.apply(lambda x: x.Uniprot_ID in uniprot_ids, axis=1)]

#     # Get metadata 
#     uni_id_lst = []
#     f_lst = []
#     gene_lst = []
#     hgnc_lst = []
#     refseq_lst = []
#     mane_missing_lst = []
#     seq_lst = []

#     list_prot_path = get_pdb_path_list_from_dir(path_to_pdb)
#     for path_structure in tqdm(list_prot_path, total=len(list_prot_path), desc="Generating sequence df"):
#         uniprot_id, f = get_af_id_from_pdb(path_structure).split("-F")
        
#         if uniprot_id in id_mapping.Uniprot_ID.unique() and uniprot_id not in uni_id_lst:
#             row = id_mapping[id_mapping["Uniprot_ID"] == uniprot_id].values[0]
#             refseq, _, mane_missing, hugo_id, hgnc_id = row
#             seq = "".join(list(get_seq_from_pdb(path_structure)))
#             uni_id_lst.append(uniprot_id)
#             f_lst.append(f)
#             gene_lst.append(hugo_id)
#             hgnc_lst.append(hgnc_id)
#             refseq_lst.append(refseq)
#             mane_missing_lst.append(mane_missing)
#             seq_lst.append(seq)

#     seq_df = pd.DataFrame({"Uniprot_ID" : uni_id_lst, 
#                            "F" : f_lst, 
#                            "Gene" : gene_lst, 
#                            "HGNC_ID" : hgnc_lst, 
#                            "Refseq_prot" : refseq_lst,
#                            "Mane_missing" : mane_missing_lst,
#                            "Seq" : seq_lst
#                            }).sort_values(["Gene", "F"])
    
#     return seq_df


#==============================================
# Get DNA sequence using EMBL backtranseq API 
# (not 100% reliable but available fo most seq)
#==============================================

def backtranseq(protein_seqs, organism = "Homo sapiens"):
    """
    Perform backtranslation from proteins to DNA sequences using EMBOS backtranseq.    
    """
    
    # Define the API endpoints
    run_url = "https://www.ebi.ac.uk/Tools/services/rest/emboss_backtranseq/run"
    status_url = "https://www.ebi.ac.uk/Tools/services/rest/emboss_backtranseq/status/"
    result_url = "https://www.ebi.ac.uk/Tools/services/rest/emboss_backtranseq/result/"
    
    # Define the parameters for the API request (an email address must be included)
    params = {"email": "stefano.pellegrini@irbbarcelona.org",
              "sequence": protein_seqs,
              "outseqformat": "plain",
              "molecule": "dna",
              "organism": organism}
    
    # Submit the job request and retrieve the job ID
    response = "INIT"
    while str(response) != "<Response [200]>":
        if response != "INIT":
            time.sleep(10)
        try:
            response = requests.post(run_url, data=params, timeout=160)
        except requests.exceptions.RequestException as e:                          
            response = "ERROR"                                                     
            logger.debug(f"Request failed: {e}")    
        
    job_id = response.text.strip()

    # Wait for the job to complete
    status = "INIT"
    while status != "FINISHED":
        time.sleep(20)
        try:
            result = requests.get(status_url + job_id, timeout=160)
            status = result.text.strip()                                              
        except requests.exceptions.RequestException as e:                          
            status = "ERROR"                                                     
            logger.debug(f"Request failed: {e}")                               

    # Retrieve the results of the job
    status = "INIT"
    while status != "FINISHED":
        try:
            result = requests.get(result_url + job_id + "/out", timeout=160)
            status = "FINISHED"
        except requests.exceptions.RequestException as e:
            status = "ERROR"                                                     
            logger.debug(f"Request failed: {e}")                               
        time.sleep(10)                                                       
        
    dna_seq = result.text.strip()
    
    return dna_seq


def batch_backtranseq(df, batch_size, organism = "Homo sapiens"):
    """
    Given a dataframe including protein sequences, it divides the 
    sequences into batches of a given size and run EMBOSS backtranseq 
    (https://www.ebi.ac.uk/Tools/st/emboss_backtranseq/) to translate 
    them into DNA sequences.
    """
    
    batches = df.groupby(df.index // batch_size)
    lst_batches = []

    # Iterate over batches
    for i, batch in tqdm(batches, total=len(batches), desc="Backtranseq"):
        
        # Avoid sending too many request
        if i+1 == 30:
            logger.debug("Reached maximum number of requests: waiting 180s..")
            time.sleep(180)
            
        # Get input format for backtranseq 
        batch_seq = "\n".join(batch.reset_index(drop=True).apply(lambda x: f'>\n{x["Seq"]}', axis=1).values)

        # Run backtranseq
        batch_dna = backtranseq(batch_seq, organism = organism)

        # Parse output
        batch_dna = re.split(">EMBOSS_\d+", batch_dna.replace("\n", ""))[1:]

        batch["Seq_dna"] = batch_dna
        lst_batches.append(batch)
        
    return pd.concat(lst_batches)


#===============================================================
# Get exons (CDS) coordinate using EMBL Proteins Coordinates API
#===============================================================

def _uniprot_request_coord(lst_uniprot_ids):
    """
    Use Coordinates from EMBL-EBI Proteins API to get 
    a json including exons coordinate and protein info.
    
    https://www.ebi.ac.uk/proteins/api/doc/#coordinatesApi
    https://doi.org/10.1093/nar/gkx237
    """
    
    prot_request = [f"{prot}" if i == 0 else f"%2C{prot}" for i, prot in enumerate(lst_uniprot_ids)]
    requestURL = f"https://www.ebi.ac.uk/proteins/api/coordinates?offset=0&size=100&accession={''.join(prot_request)}"
    
    status = "INIT"
    while status != "FINISHED":
        if status != "INIT":
            time.sleep(10)
        try:
            r = requests.get(requestURL, headers={ "Accept" : "application/json"}, timeout=160)
            if r.ok:
                status = "FINISHED"
            else:
                logger.debug(f"Error occurred after successfully sending request. Status: {r.raise_for_status()}")             
                status = "ERROR"
        except requests.exceptions.RequestException as e:                          
            status = "ERROR"                                                     
            logger.debug(f"Request failed: {e}")                               
    
    for dictio in json.loads(r.text):

        yield dictio


def get_batch_exons_coord(batch_ids):
    """
    Parse the json obtained from the Coordinates service extracting 
    exons coordinates and protein info.
    
    https://www.ebi.ac.uk/proteins/api/doc/#coordinatesApi
    https://doi.org/10.1093/nar/gkx237
    """
    
    lst_uni_id = []
    lst_ens_gene_id = []
    lst_ens_transcr_id = []
    lst_seq = []
    lst_chr = []
    lst_reverse = []
    lst_ranges = []
    
    for dic in _uniprot_request_coord(batch_ids):
        
        uni_id = dic["accession"]
        seq = dic["sequence"]
        
        # Iterate throug the transcripts
        for i in range(len(dic["gnCoordinate"])):

            ens_gene_id = dic["gnCoordinate"][i]["ensemblGeneId"]
            ens_transcr_id = dic["gnCoordinate"][i]["ensemblTranscriptId"]
            dic_loc = dic["gnCoordinate"][i]["genomicLocation"]
            exons = dic_loc["exon"]
            chromosome = dic_loc["chromosome"]
            reverse_strand = dic_loc["reverseStrand"]
            ranges = []

            for exon in exons:
                exon = exon["genomeLocation"]
                if "begin" in exon.keys():
                    begin = exon["begin"]["position"]
                    end = exon["end"]["position"]
                else:
                    begin = exon["position"]["position"]
                    end = exon["position"]["position"]
                ranges.append((begin, end))

            # Check if the DNA seq of the transcript is equal the codon lenght
            start = 0 if not int(reverse_strand) else 1
            end = 1 if not int(reverse_strand) else 0
            len_dna = sum([coord[end] - coord[start] + 1 for coord in ranges])
            if len_dna / 3 == len(seq):
                break
            # If there is no transcript with matching sequence, return NA
            elif i == len(dic["gnCoordinate"]) - 1:
                ens_gene_id = np.nan
                ens_transcr_id = np.nan
                chromosome = np.nan
                reverse_strand = np.nan
                ranges = np.nan
        
        lst_uni_id.append(uni_id)
        lst_ens_gene_id.append(ens_gene_id)
        lst_ens_transcr_id.append(ens_transcr_id)
        lst_seq.append(seq)
        lst_chr.append(chromosome)
        reverse_strand = int(reverse_strand) if isinstance(reverse_strand, int) else np.nan
        ranges = str(ranges) if isinstance(ranges, list) else np.nan
        lst_reverse.append(reverse_strand)
        lst_ranges.append(ranges)

    return pd.DataFrame({"Uniprot_ID" : lst_uni_id, 
                         "Ens_Gene_ID" : lst_ens_gene_id, 
                         "Ens_Transcr_ID" : lst_ens_transcr_id, 
                         "Seq" : lst_seq, 
                         "Chr" : lst_chr, 
                         "Reverse_strand" : lst_reverse, 
                         "Exons_coord" : lst_ranges})


def get_exons_coord(ids, batch_size=100):
    """
    Use the Coordinates service from Proteins API of EMBL-EBI to get 
    exons coordinate and proteins info of all provided Uniprot IDs.
    
    https://www.ebi.ac.uk/proteins/api/doc/#coordinatesApi
    https://doi.org/10.1093/nar/gkx237
    """
    
    lst_df = []
    batches_ids = [ids[i:i+batch_size] for i in range(0, len(ids), batch_size)]

    for batch_ids in tqdm(batches_ids, total=len(batches_ids), desc="Adding exons coordinate"):
        
        batch_df = get_batch_exons_coord(batch_ids)
        
        # Identify unmapped IDs and add them as NaN rows
        unmapped_ids = list(set(batch_ids).difference(set(batch_df.Uniprot_ID.unique())))
        nan = np.repeat(np.nan, len(unmapped_ids))
        nan_rows = pd.DataFrame({'Uniprot_ID' : unmapped_ids,
                                 'Ens_Gene_ID' : nan,
                                 'Ens_Transcr_ID' : nan,
                                 'Seq': nan,
                                 'Chr': nan,
                                 'Reverse_strand' : nan,
                                 'Exons_coord': nan})

        batch_df = pd.concat([batch_df, nan_rows], ignore_index=True)
        lst_df.append(batch_df)
    
    return pd.concat(lst_df).reset_index(drop=True)


#=======================================================================
# Get DNA sequence and trin-context of reference genome from coordinates
#=======================================================================

def per_site_trinucleotide_context(seq, no_flanks=False):
    """
    Given a DNA sequence rapresenting an exon (CDS) with -1 and +1 flanking 
    sites, return a list including the trinucletide context of each site.
    """
    
    # If there is no info about flanking regions, add most the two most common ones
    if no_flanks:
        seq = f"C{seq}T"
    
    return [f"{seq[i-1]}{seq[i]}{seq[i+1]}" for i in range(1, len(seq) - 1)]


def get_ref_dna_and_context(row, genome_fun):
    """
    Given a row of the sequence dataframe including exons (CDS) coordinate, 
    get the corresponding DNA sequence and the trinucleotide context of each 
    site of the sequence taking into account flaking regions at splicing sites.
    """
    
    transcript_id = row["Ens_Transcr_ID"]
    
    try:
        lst_exon_coord = ast.literal_eval(row["Exons_coord"])
        reverse = row["Reverse_strand"]
        chrom = row["Chr"]

        lst_exon_tri_context = []
        lst_exon_seq = []
        
        for region in lst_exon_coord:
            
            # Retrieve reference seq of the exon (CDS)
            start = 0 if not reverse else 1
            end = 1 if not reverse else 0
            segment_len = region[end] - region[start] + 1
            seq_with_flanks = genome_fun(chrom, region[start] - 1, size = segment_len + 2)
            
            # Get reverse complement if reverse strand
            if reverse:
                seq_with_flanks = str(Seq(seq_with_flanks).reverse_complement())
            
            # Get the trinucleotide context of each site of the exon (CDS)
            per_site_tri_context = ",".join(per_site_trinucleotide_context(seq_with_flanks))
            
            # Get the reference DNA sequence
            seq = seq_with_flanks[1:-1]

            lst_exon_tri_context.append(per_site_tri_context)
            lst_exon_seq.append(seq)

        return transcript_id, "".join(lst_exon_seq), ",".join(lst_exon_tri_context), 1
    
    except Exception as e:
        if not str(e).startswith("Sequence"):
            logger.warning(f"Error occurred during retrieving DNA seq from ref coordinates {transcript_id}: {e}")
        
        return transcript_id, np.nan, np.nan, -1


def add_ref_dna_and_context(seq_df, genome_fun=hg38):
    """
    Given as input the sequence dataframe including exons (CDS) coordinate, 
    add the corresponding DNA sequence and the trinucleotide context of each 
    site of the sequence taking into account flaking regions at splicing sites.
    """
    
    seq_df = seq_df.copy()
    seq_df_with_coord = seq_df[["Ens_Transcr_ID", "Chr", "Reverse_strand", "Exons_coord"]].dropna().drop_duplicates()
    ref_dna_and_context = seq_df_with_coord.apply(lambda x: get_ref_dna_and_context(x, genome_fun), axis=1, result_type='expand')
    ref_dna_and_context.columns = ["Ens_Transcr_ID", "Seq_dna", "Tri_context", "Reference_info"]
    seq_df = seq_df.merge(ref_dna_and_context, on=["Ens_Transcr_ID"], how="left")
    seq_df["Reference_info"] = seq_df["Reference_info"].fillna(-1).astype(int)
    
    return seq_df


#======
# Utils
#======

# def asssess_similarity(seq_df, on="all"):
#     """
#     Assess similarity between protein sequences and translated 
#     DNA sequences by reference coordinates or/and backtranslation.
#     """
    
#     logger.debug("Assessing similarity..")
    
#     # Ref seq and backtranseq backtranslation
#     if on == "all":
#         seq_df["Seq_similarity"] = seq_df.apply(lambda x: get_seq_similarity(x.Seq, translate_dna(x.Seq_dna)), axis=1) 
#         avg_sim = np.mean(seq_df["Seq_similarity"].dropna())
#         seq_ref = seq_df[seq_df["Reference_info"] == 1]
#         seq_ref = seq_ref.drop_duplicates('Uniprot_ID')
#         seq_backtr = seq_df[seq_df["Reference_info"] == 0]
#         seq_backtr = seq_backtr.drop_duplicates('Uniprot_ID')
#         avg_sim_ref = np.mean(seq_ref["Seq_similarity"].dropna())
#         avg_sim_backtr = np.mean(seq_backtr["Seq_similarity"].dropna())
#         if len(seq_ref) > 0: 
#             logger.debug(f"Average similarity ref coord backtranslation: {avg_sim_ref:.2f} for {len(seq_ref)} sequences.")      
#         else:
#             logger.debug("Ref coord backtranslated sequences not available.")
#         if len(seq_backtr) > 0: 
#             logger.debug(f"Average similarity backtranseq backtranslation: {avg_sim_backtr:.2f} for {len(seq_backtr)} sequences.") 
#         else:
#             logger.debug("Backtranseq backtranslated sequences not available.")
#         if avg_sim < 1:                       
#             logger.warning(f"Mismatch between protein and translated DNA. Overall average similatity: {avg_sim:.2f}.")
#             logger.warning("3D-clustering analysis using mutational profile and/or mutability can be severely affected.")
#         else:
#             logger.debug("Final similarity check: PASS")
        
#         return seq_df
         
#     # Ref seq backtranslation only
#     else:
#         seq_df["Seq_similarity"] = seq_df.apply(lambda x: get_seq_similarity(x.Seq, translate_dna(x.Seq_dna)), axis=1)
#         avg_sim = np.mean(seq_df["Seq_similarity"].dropna())
#         if avg_sim < 1:
#             tot_mismatch = sum(seq_df['Seq_similarity'].dropna() < 1)
#             tot_seq = len(seq_df['Seq_similarity'].dropna())
#             mismatch_ratio = tot_mismatch / tot_seq
#             logger.warning(f"Mismatch between proteins and translated DNA by reference coordinates on {tot_mismatch} of {tot_seq}.")
#             logger.warning(f"Protein and translated DNA mean similarity < 1: {avg_sim:.2f}")
#             if mismatch_ratio > 0.4:
#                 logger.warning(f"Mismatch ratio > 0.4: {mismatch_ratio:.2f}. Dropping all {tot_mismatch} ref DNA sequnces..") 
#                 seq_df["Reference_info"] = 0
#             else:
#                 logger.warning(f"Dropping {tot_mismatch} mismatching ref DNA sequnces..") 
#                 seq_df.loc[seq_df['Seq_similarity'] < 1, "Reference_info"] = 0
#         else:
#             logger.debug("Ref seq similarity check: PASS")
                
#         return seq_df
                

# def drop_gene_duplicates(df):
#     """
#     Issue: multiple Uniprot IDs might map to the same HUGO symbol.
#     Drop Uniprot ID of gene duplicates by prioritizing the following order:
#     - Uniprot ID whose sequence has been obtained using reference info
#     - Uniprot ID having perfect match with MANE transcripts
#     """

#     df = df.copy()
#     df_ref_1 = df[df['Reference_info'] == 1]
#     df_ref_0 = df[df['Reference_info'] == 0]
    
#     # Prioritize rows where MANE not NaN for genes with Reference_info (or just keep first one)
#     df_ref_1 = df_ref_1.sort_values(by='Refseq_prot', ascending=False).drop_duplicates(subset='Gene')
    
#     # Same as step 1 but for genes without Reference_info
#     df_ref_0 = df_ref_0.sort_values(by='Refseq_prot', ascending=False).drop_duplicates(subset='Gene')
#     df_ref_0 = df_ref_0.drop_duplicates(subset='Gene', keep='first')
    
#     # Concat and prioritize rows with Reference_info (or just keep first one)
#     result_df = pd.concat([df_ref_1, df_ref_0]).sort_values(["Gene", "Reference_info"], ascending=False)
#     result_df = result_df.drop_duplicates(subset='Gene', keep='first').sort_values(["Gene"]) 
#     result_df = result_df.reset_index(drop=True)
    
#     # Check if there are still duplicates
#     n_duplicates = sum(result_df["Gene"].value_counts() > 1)
#     if n_duplicates > 0:
#         logger.warning(f"Found {n_duplicates} duplicates gene entries: Mapping HUGO Symbol to protein info might be affected.")
#     else:
#         logger.debug("Duplicates gene entries correctly removed!")

#     return result_df.reset_index(drop=True)

                
#=========
# WRAPPERS
#=========

# TODO: test if with the added "if gene not in" is avoiding duplicates Hugo Symbol with different Uni IDs
#       - In general test if there are multiple rows with same gene name in the final seq_df
#       - It seems that backtranseq is used for Uniprot ID connected to genes that were already processed with ref DNA: must avoid it

def add_extra_genes_to_seq_df(seq_df, uniprot_to_gene_dict):
    """
    If multiple genes are mapping to a given Uniprot_ID, add 
    each gene name with corresponding sequence info to the seq_df.
    """
    
    lst_added_genes = []
    lst_extra_genes_rows = []
    for _, seq_row in seq_df.iterrows():
        uni_id = seq_row["Uniprot_ID"]
        if uni_id in uniprot_to_gene_dict:
            gene_id = uniprot_to_gene_dict[uni_id]
            
            if not pd.isna(gene_id):
                gene_id = gene_id.split()
                
                if len(gene_id) > 1:
                    for gene in gene_id:
                        if gene != seq_row["Gene"] and gene not in lst_added_genes:
                            row = seq_row.copy()
                            row["Gene"] = gene
                            lst_extra_genes_rows.append(row)
                            lst_added_genes.append(gene)
    
    seq_df_extra_genes = pd.concat(lst_extra_genes_rows, axis=1).T
    
    # Remove rows with multiple symbols and drop duplicated ones
    seq_df = pd.concat((seq_df, seq_df_extra_genes))
    seq_df = seq_df.dropna(subset=["Gene"])
    seq_df = seq_df[seq_df.apply(lambda x: len(x["Gene"].split()), axis =1) == 1].reset_index(drop=True)
    # seq_df.Uniprot_ID = seq_df.Uniprot_ID.str.replace(";", "")
    # seq_df.Gene = seq_df.Gene.str.replace(";", "")
    seq_df = seq_df.drop_duplicates().reset_index(drop=True)
    
    return seq_df


# def get_seq_df(datasets_dir, 
#                output_seq_df, 
#                uniprot_to_gene_dict = None, 
#                organism = "Homo sapiens",
#                mane=False):
#     """
#     Generate a dataframe including IDs mapping information (Gene and Uniprot_ID),
#     AlphaFold 2 fragment number (F); the protein (Seq) and DNA sequence (DNA_seq) 
#     obtained by EMBL backtranseq, as well as the genomic coordinate of the exons 
#     (Chr and Exons_coord), which are used to compute the per-residue probability 
#     of missense mutation. Also, for a subset of genes, it include the DNA sequence 
#     of the reference genome (Seq_dna_ref) and its per-site trinucleotide context
#     taing into account flanking regions at splicing sites (Tri_context).
#     """

#     # if mane == False:
#     #     # Load Uniprot ID to HUGO mapping
#     #     pdb_dir = os.path.join(datasets_dir, "pdb_structures")
#     #     uniprot_ids = os.listdir(pdb_dir)
#     #     uniprot_ids = [uni_id.split("-")[1] for uni_id in list(set(uniprot_ids)) if ".pdb" in uni_id]
#     #     if uniprot_to_gene_dict is not None and uniprot_to_gene_dict != "None":
#     #         uniprot_to_gene_dict = json.load(open(uniprot_to_gene_dict)) 
#     #     else:
#     #         logger.debug("Retrieving Uniprot_ID to Hugo symbol mapping information..")
#     #         uniprot_to_gene_dict = uniprot_to_hugo(uniprot_ids)  
    
#     # Load Uniprot ID to HUGO mapping
#     pdb_dir = os.path.join(datasets_dir, "pdb_structures")                                                                  ###### -X
#     uniprot_ids = os.listdir(pdb_dir)                                                                                       ###### -X
#     uniprot_ids = [uni_id.split("-")[1] for uni_id in list(set(uniprot_ids)) if ".pdb" in uni_id]                           ###### -X
#     logger.debug("Retrieving Uniprot_ID to Hugo symbol mapping information..")                                              ###### -X
#     uniprot_to_gene_dict = uniprot_to_hugo(uniprot_ids)                                                                     ###### -X

#     # # Create a dataframe with protein sequences
#     # logger.debug("Generating sequence df..")
#     # if mane == False:  
#     #     seq_df = initialize_seq_df(pdb_dir, uniprot_to_gene_dict)
#     #     seq_df["Refseq_prot"] = np.nan
#     # else:
#     #     seq_df = initialize_seq_df_from_mane(datasets_dir) 
#     logger.debug("Generating sequence df..")                                                                                ###### -X
#     seq_df = initialize_seq_df(pdb_dir, uniprot_to_gene_dict)                                                              ###### -X
    
#     # Add coordinates for mutability integration (entries in Proteins API)
#     logger.debug("Retrieving exons coordinate..")
#     coord_df = get_exons_coord(seq_df["Uniprot_ID"].unique())
#     seq_df = seq_df.merge(coord_df, on=["Seq", "Uniprot_ID"], how="left").reset_index(drop=True)
    
#     # Add ref DNA seq and its per-site trinucleotide context (entries in Proteins API)
#     logger.debug("Retrieving exons DNA sequence from reference genome..")
#     if organism == "Homo sapiens":
#         genome_fun = hg38
#     elif organism == "Mus musculus":
#         genome_fun = mm10
#     else:
#         raise RuntimeError(f"Failed to recognize '{organism}' as organism. Currently accepted ones are 'Homo sapiens' and 'Mus musculus'. Exiting..")
#     seq_df = add_ref_dna_and_context(seq_df, genome_fun)
    
#     # # Check if ref DNA sequences match protein sequences
#     # seq_df = asssess_similarity(seq_df, on="ref")
    
#     # Add DNA seq for genes with available MANE-associated structure
#     # TO DO: load mane mapping, retrieve transcript ID, use Ensembl sequence API to get sequence of DNA
    
#     # Add DNA seq for any other genes 
#     logger.debug("Performing backtranseq backtranslation for genes without available ref DNA..")
#     seq_df_backtranseq = seq_df[seq_df["Reference_info"] == 0].reset_index(drop=True)    
#     seq_df_backtranseq = batch_backtranseq(seq_df_backtranseq, 500, organism=organism)  
    
#     # Get trinucleotide context for (entries not in Proteins API)
#     seq_df_backtranseq["Tri_context"] = seq_df_backtranseq["Seq_dna"].apply(
#         lambda x: ",".join(per_site_trinucleotide_context(x, no_flanks=True)))
#     seq_df = pd.concat((seq_df[seq_df["Reference_info"] == 1], seq_df_backtranseq)
#                        ).sort_values("Uniprot_ID").reset_index(drop=True)
    
#     # # Add multiple genes mapping to the same Uniprot_ID
#     # if mane == False:
#     #     seq_df = add_extra_genes_to_seq_df(seq_df, uniprot_to_gene_dict)
    
#     # # # Assess overal similarity
#     # # seq_df = asssess_similarity(seq_df, on="all")
        
#     # # Drop duplicates
#     # seq_df = drop_gene_duplicates(seq_df)
        
#     # Save
#     seq_df.to_csv(output_seq_df, index=False, sep="\t")                    
#     logger.debug(f"Sequences dataframe saved in: {output_seq_df}")


def get_mane_summary(path_to_file, v=1.0, max_attempts=15):
    
    mane_summary_url = f"https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_{v}/MANE.GRCh38.v{v}.summary.txt.gz"
    logger.debug(f"Downloading MANE summary file {mane_summary_url} to {path_to_file}")
    attempts = 0
    
    while not os.path.exists(path_to_file):
        download_single_file(mane_summary_url, path_to_file, threads=1)
        attempts += 1
        if attempts >= max_attempts:
            raise RuntimeError(f"Failed to download MANE summary file after {max_attempts} attempts. Exiting...")
        time.sleep(5)


def select_uni_id(ids_tuple, all_ids):
    """
    Return the first Uniprot ID present in the list of IDs 
    (list of structures available). If no ID of the tuple 
    maps a downloaded structure, return NA.
    """

    for uni_id in ids_tuple:
        if uni_id in all_ids:
            return uni_id

    return np.nan

    
def get_mane_to_af_mapping(datasets_dir, include_not_af=False, mane_version=1.0):

    mane_to_af = pd.read_csv(os.path.join(datasets_dir, "mane_refseq_prot_to_alphafold.csv"))
    mane_to_af = mane_to_af.rename(columns={"refseq_prot" : "Refseq_prot",
                                            "uniprot_accession" : "Uniprot_ID"}).drop(columns=["alphafold"])
    path_mane_summary = os.path.join(datasets_dir, "mane_summary.txt.gz")
    if not os.path.exists(path_mane_summary):
        get_mane_summary(path_mane_summary, mane_version)
    mane = pd.read_csv(path_mane_summary, compression='gzip', sep="\t").dropna(subset=["symbol", "HGNC_ID"])
    mane = mane.rename(columns={"symbol" : "Gene",
                                "RefSeq_prot" : "Refseq_prot",
                                "Ensembl_Gene" : "Ens_Gene_ID",
                                "Ensembl_nuc" : "Ens_Transcr_ID",
                                "GRCh38_chr" : "Chr",
                                "chr_strand" : "Reverse_strand"})
    #mane = mane[["Gene", "HGNC_ID", "Refseq_prot", "Ens_Gene_ID", "Ens_Transcr_ID", "Chr", "Reverse_strand"]]
    mane = mane[["Gene", "Refseq_prot", "Ens_Gene_ID", "Ens_Transcr_ID", "Chr", "Reverse_strand"]]
    mane_mapping = mane_to_af.merge(mane, how="left", on="Refseq_prot").dropna()
    mane_mapping.Reverse_strand = mane_mapping.Reverse_strand.map({"+" : 0, "-" : 1})
    mane_mapping.Ens_Gene_ID = mane_mapping.Ens_Gene_ID.apply(lambda x: x.split(".")[0])
    mane_mapping.Ens_Transcr_ID = mane_mapping.Ens_Transcr_ID.apply(lambda x: x.split(".")[0])

    # Select first Uniprot ID if multiple ones are present
    mane_mapping["Uniprot_ID"] = mane_mapping.apply(lambda x: 
                                            select_uni_id(x.Uniprot_ID.split(";"), uniprot_ids) if len(x.Uniprot_ID.split(";")) > 1 
                                            else x.Uniprot_ID, axis=1)
    mane_mapping = mane_mapping.reset_index(drop=True)
    
    if include_not_af:
        mane_not_af = mane[~mane.Gene.isin(mane_mapping.Gene)].reset_index(drop=True)
        
        return mane_mapping, mane_not_af
        
    else:
        return mane_mapping


def get_ref_dna_from_ensembl(transcript_id):
    """
    Use Ensembl GET sequence rest API to obtain CDS DNA 
    sequence from Ensembl transcript ID.
    
    https://rest.ensembl.org/documentation/info/sequence_id
    """

    server = "https://rest.ensembl.org"
    ext = f"/sequence/id/{transcript_id}?type=cds"
    
    status = "INIT"
    i = 0
    while status != "FINISHED":

        try:
            r = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"}, timeout=160)
            if not r.ok:
                r.raise_for_status()
                status = "ERROR"
                sys.exit()
            else:
                status = "FINISHED"
            
        except requests.exceptions.RequestException as e:
            i += 1
            status = "ERROR"
            if i == 100:
                logger.debug(f"Failed to retrieve sequence for {transcript_id} ({e}): Skipping..")
                return np.nan
                
            time.sleep(0.5)

    seq_dna = "".join(r.text.strip().split("\n")[1:])

    return seq_dna[:len(seq_dna)-3]


def get_ref_dna_from_ensembl_wrapper(ensembl_id):
    """
    Wrapper for multiple processing function using 
    Ensembl Get sequence rest API.
    """
    
    return get_ref_dna_from_ensembl(ensembl_id)


def get_ref_dna_from_ensembl_mp(seq_df, cores):
    """
    Multiple processing function to use Ensembl GET sequence 
    rest API to obtain CDS DNA sequence from Ensembl transcript ID.
    
    https://rest.ensembl.org/documentation/info/sequence_id
    """

    pool = multiprocessing.Pool(processes=cores)
    seq_df = seq_df.copy()
    seq_df["Seq_dna"] = pool.map(get_ref_dna_from_ensembl_wrapper, seq_df.Ens_Transcr_ID)
    pool.close()
    pool.join()

    return seq_df


def drop_gene_duplicates(df):
    """
    Issue: multiple Uniprot IDs might map to the same HUGO symbol.
    Drop Uniprot ID of gene duplicates by prioritizing reference 
    info status (Proteins API > Ensembl API > Backtranseq API).
    """
    
    df = df.copy()
    df = df.sort_values(["Gene", "Reference_info"], ascending=False).drop_duplicates(subset='Gene').reset_index(drop=True)
    df.Uniprot_ID = df.Uniprot_ID.str.replace(";", "")
    df.Gene = df.Gene.str.replace(";", "")

    # Check if there are still duplicates
    n_duplicates = sum(df["Gene"].value_counts() > 1)
    if n_duplicates > 0:
        logger.debug(f"Found {n_duplicates} duplicates gene entries: Mapping HUGO Symbol to protein info might be affected.")
    else:
        logger.debug("Duplicates gene entries correctly removed!")

    return df.reset_index(drop=True)
    

def process_seq_df(seq_df, datasets_dir, organism, uniprot_to_gene_dict, num_cores=1, rm_weird_chr=False):
    """
    Retrieve DNA sequence and tri-nucleotide context 
    for each structure in the initialized dataframe 
    prioritizing structures obtained from transcripts 
    whose exon coordinates are available in the Proteins API.
    
    Reference_info labels:
        1 : Transcript ID, exons coord, seq DNA obtained from Proteins API
        0 : Transcript ID retrieved from MANE and seq DNA from Ensembl 
       -1 : Not available transcripts, seq DNA retrieved from Backtranseq API
    """
    
    # Process entries in Proteins API (Reference_info 1)
    #---------------------------------------------------
    
    # Add coordinates for mutability integration (entries in Proteins API)
    logger.debug(f"Retrieving CDS DNA seq from reference genome (Proteins API): {len(seq_df['Uniprot_ID'].unique())} structures..")
    coord_df = get_exons_coord(seq_df["Uniprot_ID"].unique())
    seq_df = seq_df.merge(coord_df, on=["Seq", "Uniprot_ID"], how="left").reset_index(drop=True)
    
    # Add ref DNA seq and its per-site trinucleotide context (entries in Proteins API)
    if organism == "Homo sapiens":
        genome_fun = hg38
    elif organism == "Mus musculus":
        genome_fun = mm10
    else:
        raise RuntimeError(f"Failed to recognize '{organism}' as organism. Currently accepted ones are 'Homo sapiens' and 'Mus musculus'. Exiting..")
    seq_df = add_ref_dna_and_context(seq_df, genome_fun)
    seq_df_tr = seq_df[seq_df["Reference_info"] == 1]
    seq_df_notr = seq_df[seq_df["Reference_info"] == -1]
    
    
    # Process entries not in Proteins API (Reference_info 0 & -1)
    #------------------------------------------------------------
    
    # Retrieve transcripts from MANE metadata
    if organism == "Homo sapiens":                                                 
        mane_mapping = get_mane_to_af_mapping(datasets_dir)
        seq_df_notr = seq_df_notr.drop(columns=["Ens_Gene_ID", "Ens_Transcr_ID", "Chr", "Reverse_strand"])
        seq_df_notr = seq_df_notr.merge(mane_mapping.drop(columns=["Gene", "Refseq_prot"]), on="Uniprot_ID", how="left").dropna(subset="Gene")

        # Assign reference label to transcripts retrieved from MANE
        seq_df_natr = seq_df_notr[seq_df_notr.Ens_Transcr_ID.isna()]
        seq_df_manetr = seq_df_notr.copy()
        seq_df_manetr = seq_df_manetr.dropna(subset="Ens_Transcr_ID")
        seq_df_manetr["Reference_info"] = 0
        
        # Remove genes from weird chrs
        if rm_weird_chr:
            chrs_lst = [f"{i}" for i in range(1, 23)] + ['X', 'Y']
            seq_df_manetr = seq_df_manetr[seq_df_manetr['Chr'].isna() | seq_df_manetr['Chr'].isin(chrs_lst)]

        # Add DNA seq from Ensembl for structures with available transcript ID
        logger.debug(f"Retrieving CDS DNA seq from transcript ID (Ensembl API): {len(seq_df_manetr)} structures..")
        seq_df_manetr = get_ref_dna_from_ensembl_mp(seq_df_manetr, cores=num_cores)

        # Set failed and len-mismatching entries as no-transcripts entries
        failed_ix = seq_df_mane.apply(lambda x: True if pd.isna(x.Seq_dna) else len(x.Seq_dna) / 3 != len(x.Seq), axis=1)
        if sum(failed_ix) > 0:
            seq_df_manetr_failed = seq_df_manetr[failed_ix]
            seq_df_manetr = seq_df_manetr[~failed_ix]
            seq_df_manetr_failed["Reference_info"] = -1
            seq_df_manetr_failed.drop(columns=["Ens_Gene_ID", "Ens_Transcr_ID", "Reverse_strand", "Chr"])
            seq_df_notr = pd.concat((seq_df_notr, seq_df_manetr_failed))
    else:
        seq_df_manetr = pd.DataFrame()
    
    # Add DNA seq from Backtranseq for any other entry
    logger.debug(f"Retrieving CDS DNA seq for entries without available transcript ID (Backtranseq API): {len(seq_df_notr)} structures..") 
    seq_df_notr = batch_backtranseq(seq_df_notr, 500, organism=organism)  
    
    # Get trinucleotide context
    seq_df_not_uniprot = pd.concat((seq_df_manetr, seq_df_notr))
    seq_df_not_uniprot["Tri_context"] = seq_df_not_uniprot["Seq_dna"].apply(
        lambda x: ",".join(per_site_trinucleotide_context(x, no_flanks=True)))
    
    
    # Prepare final output
    #---------------------
    
    # Concat the dfs, expand multiple genes associated to the same structure, keep only one structure for each gene
    seq_df = pd.concat((seq_df_tr, seq_df_not_uniprot)).reset_index(drop=True)
    logger_report = ", ".join([f"{v}: {c}" for (v, c) in zip(seq_df.Reference_info.value_counts().index, 
                                                              seq_df.Reference_info.value_counts().values)])
    logger.debug(f"Built of sequence dataframe completed. Retrieved {len(seq_df)} structures: {logger_report}")
    seq_df.to_csv(datasets_dir + "/seq_df_pre_add_extra.tsv", index=False, sep="\t")                                                      ## TO DEL                 
    seq_df = add_extra_genes_to_seq_df(seq_df, uniprot_to_gene_dict)
    seq_df.to_csv(datasets_dir + "/seq_df_pre_drop_dupl.tsv", index=False, sep="\t")                                                      ## TO DEL               
    seq_df = drop_gene_duplicates(seq_df)
    
    return seq_df


def process_seq_df_mane(seq_df, datasets_dir, uniprot_to_gene_dict, num_cores=1):
    """
    Retrieve DNA sequence and tri-nucleotide context 
    for each structure in the initialized dataframe
    prioritizing MANE associated structures and metadata.
    
    Reference_info labels:
        1 : Transcript ID, exons coord, seq DNA obtained from Proteins API
        0 : Transcript ID retrieved from MANE and seq DNA from Ensembl 
       -1 : Not available transcripts, seq DNA retrieved from Backtranseq API
    """
    
    mane_mapping, mane_mapping_not_af = get_mane_to_af_mapping(datasets_dir, include_not_af=True)
    seq_df_mane = seq_df[seq_df.Uniprot_ID.isin(mane_mapping.Uniprot_ID)].reset_index(drop=True)
    seq_df_nomane = seq_df[~seq_df.Uniprot_ID.isin(mane_mapping.Uniprot_ID)].reset_index(drop=True)
    
    # Seq df MANE
    seq_df_mane = seq_df_mane.drop(columns=["Gene"]).merge(mane_mapping, how="left", on="Uniprot_ID")
    seq_df_mane["Reference_info"] = 0
    
    # Add DNA seq from Ensembl for structures with available transcript ID
    logger.debug(f"Retrieving CDS DNA seq from transcript ID (Ensembl API): {len(seq_df_mane)} structures..")
    seq_df_mane = get_ref_dna_from_ensembl_mp(seq_df_mane, cores=num_cores)
    
    # Set failed and len-mismatching entries as no-transcripts entries
    failed_ix = seq_df_mane.apply(lambda x: True if pd.isna(x.Seq_dna) else len(x.Seq_dna) / 3 != len(x.Seq), axis=1)
    if sum(failed_ix) > 0:
        seq_df_mane_failed = seq_df_mane[failed_ix]
        seq_df_mane = seq_df_mane[~failed_ix]
        seq_df_mane_failed = seq_df_mane_failed.drop(columns=["Ens_Gene_ID", "Ens_Transcr_ID", "Reverse_strand", 
                                                              "Chr", "Refseq_prot", "Reference_info", "Seq_dna"])
        seq_df_nomane = pd.concat((seq_df_nomane, seq_df_mane_failed))


    # Seq df not MANE
    seq_df_nomane = add_extra_genes_to_seq_df(seq_df_nomane, uniprot_to_gene_dict)
    seq_df_nomane = seq_df_nomane[seq_df_nomane.Gene.isin(mane_mapping_not_af.Gene)]
    
    # Retrieve seq from coordinates
    logger.debug(f"Retrieving CDS DNA seq from reference genome (Proteins API): {len(seq_df_nomane['Uniprot_ID'].unique())} structures..")
    coord_df = get_exons_coord(seq_df_nomane["Uniprot_ID"].unique())
    seq_df_nomane = seq_df_nomane.merge(coord_df, on=["Seq", "Uniprot_ID"], how="left").reset_index(drop=True)
    seq_df_nomane = add_ref_dna_and_context(seq_df_nomane, hg38)                                                        ## GOT ERROR HERE
    seq_df_nomane_tr = seq_df_nomane[seq_df_nomane["Reference_info"] == 1]
    seq_df_nomane_notr = seq_df_nomane[seq_df_nomane["Reference_info"] == -1]
    
    # Add DNA seq from Backtranseq for any other entry
    logger.debug(f"Retrieving CDS DNA seq for genes without available transcript ID (Backtranseq API): {len(seq_df_nomane_notr)} structures..") 
    seq_df_nomane_notr = batch_backtranseq(seq_df_nomane_notr, 500, organism="Homo sapiens")  
    
    # Get trinucleotide context
    seq_df_not_uniprot = pd.concat((seq_df_mane, seq_df_nomane_notr))
    seq_df_not_uniprot["Tri_context"] = seq_df_not_uniprot["Seq_dna"].apply(
        lambda x: ",".join(per_site_trinucleotide_context(x, no_flanks=True)))
    
    # Prepare final output
    seq_df_nomane_tr = drop_gene_duplicates(seq_df_nomane_tr)          
    seq_df = pd.concat((seq_df_not_uniprot, seq_df_nomane_tr)).reset_index(drop=True)
    report_df = seq_df.Reference_info.value_counts().reset_index()
    report_df = report_df.rename(columns={"index" : "Source"})
    report_df.Source = report_df.Source.map({1 : "Proteins API",
                                             0 : "MANE + Ensembl",
                                            -1 : "Backtranseq"})
    logger_report = ", ".join([f"{v}: {c}" for (v, c) in zip(report_df.Source, 
                                                              report_df.Reference_info)])
    logger.debug(f"Built of sequence dataframe completed. Retrieved {len(seq_df)} structures: {logger_report}")
    
    return seq_df
    

def get_seq_df(datasets_dir, 
               output_seq_df, 
               uniprot_to_gene_dict = None, 
               organism = "Homo sapiens",
               mane=False,
               num_cores=1,
               rm_weird_chr=False):
    """
    Generate a dataframe including IDs mapping information, the protein 
    sequence, the DNA sequence and its tri-nucleotide context, which is 
    used to compute the per-residue probability of missense mutation. 
    
    The DNA sequence and the tri-nucleotide context can be obtained by the
    Proteins API (obtaining the coordinates of the exons and then using them 
    to obtain the sequence from the reference genome), Ensembl GET sequence API
    (from Ensembl transcript ID), and from Backtranseq API for all other entries.
    
    https://www.ebi.ac.uk/proteins/api/doc/
    https://rest.ensembl.org/documentation/info/sequence_id
    https://www.ebi.ac.uk/jdispatcher/st/emboss_backtranseq
    """
    
    # Initialization
    #===============
    
    # Load Uniprot ID to HUGO and MANE to AF mapping
    pdb_dir = os.path.join(datasets_dir, "pdb_structures")                                                                
    uniprot_ids = os.listdir(pdb_dir)                                                                                  
    uniprot_ids = [uni_id.split("-")[1] for uni_id in list(set(uniprot_ids)) if ".pdb" in uni_id]                      
    logger.debug("Retrieving Uniprot ID to HUGO symbol mapping information..")                                         
    uniprot_to_gene_dict = uniprot_to_hugo(uniprot_ids)       
    # Workaround if the direct request to UniprotKB stops working (it has happened temporarily)
    if all(pd.isna(k) for k in uniprot_to_gene_dict.keys()):
        logger.warning(f"Failed to retrieve Uniprot ID to HUGO symbol mapping directly from UniprotKB.")
        logger.warning(f"Retrying using Unipressed API client (only first HUGO symbol entry will be mapped)..")
        uniprot_to_gene_dict = uniprot_to_hugo(uniprot_ids)   
    
    # Create a dataframe with protein sequences
    logger.debug("Initializing sequence df..")                                                                            
    seq_df = initialize_seq_df(pdb_dir, uniprot_to_gene_dict)                                                            
    
    if mane:
       seq_df = process_seq_df_mane(seq_df, datasets_dir, uniprot_to_gene_dict, num_cores)
    else:
        seq_df = process_seq_df(seq_df, datasets_dir, organism, uniprot_to_gene_dict, num_cores, rm_weird_chr)
    
    # Save
    seq_df.to_csv(output_seq_df, index=False, sep="\t")                    
    logger.debug(f"Sequences dataframe saved in: {output_seq_df}")