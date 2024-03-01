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
import requests
import numpy as np
import daiquiri
import pandas as pd
from tqdm import tqdm
from bgreference import hg38, mm10
from Bio.Seq import Seq

from scripts import __logger_name__
from scripts.datasets.utils import (get_af_id_from_pdb,
                                    get_pdb_path_list_from_dir,
                                    get_seq_from_pdb, get_seq_similarity,
                                    translate_dna, uniprot_to_hugo)

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
    
    # Get all PDB path in directory whose IDs can be mapped to HUGO symbols
    list_prot_path = get_pdb_path_list_from_dir(input_path)
    list_prot_path = [path for path in list_prot_path if get_af_id_from_pdb(path).split("-F")[0] in uniprot_to_gene_dict.keys()]
    pdb_not_mapped = set([get_af_id_from_pdb(path).split("-F")[0] 
                          for path in list_prot_path if get_af_id_from_pdb(path).split("-F")[0] not in uniprot_to_gene_dict.keys()])
    if len(pdb_not_mapped) > 0:                                 
        logger.warning(f"{len(pdb_not_mapped)} Uniprot-ID not found in the Uniprot-HUGO mapping dictionary.")

    # Get Uniprot ID, HUGO, F and protein sequence of any PDB in dir
    gene_lst = []
    uni_id_lst = []
    f_lst = []
    seq_lst = []

    for path_structure in tqdm(list_prot_path, total=len(list_prot_path), desc="Generating sequence df"):
        uniprot_id, f = get_af_id_from_pdb(path_structure).split("-F")
        gene = uniprot_to_gene_dict[uniprot_id]
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
            logger.debug("Reached maximum number of requests: waiting 180s...")
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
    lst_ref_mismatch = []
    
    for dic in _uniprot_request_coord(batch_ids):
        
        uni_id = dic["accession"]
        seq = dic["sequence"]
        ref_mismatch = 0
        
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
                ref_mismatch = 1
        
        lst_uni_id.append(uni_id)
        lst_ens_gene_id.append(ens_gene_id)
        lst_ens_transcr_id.append(ens_transcr_id)
        lst_seq.append(seq)
        lst_chr.append(chromosome)
        reverse_strand = int(reverse_strand) if isinstance(reverse_strand, int) else np.nan
        ranges = str(ranges) if isinstance(ranges, list) else np.nan
        lst_reverse.append(reverse_strand)
        lst_ranges.append(ranges)
        lst_ref_mismatch.append(ref_mismatch)

    return pd.DataFrame({"Uniprot_ID" : lst_uni_id, 
                         "Ens_Gene_ID" : lst_ens_gene_id, 
                         "Ens_Transcr_ID" : lst_ens_transcr_id, 
                         "Seq" : lst_seq, 
                         "Chr" : lst_chr, 
                         "Reverse_strand" : lst_reverse, 
                         "Exons_coord" : lst_ranges,
                         "Ref_mismatch" : lst_ref_mismatch})


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
                                 'Exons_coord': nan,
                                 "Ref_mismatch" : nan})

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
            logger.warning(f"Error occurred during backtranslation by ref coordinates {transcript_id}: {e}")
        
        return transcript_id, np.nan, np.nan, 0


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
    seq_df["Reference_info"] = seq_df["Reference_info"].fillna(0).astype(int)
    
    return seq_df


#======
# Utils
#======

def asssess_similarity(seq_df, on="all"):
    """
    Assess similarity between protein sequences and translated 
    DNA sequences by reference coordinates or/and backtranslation.
    """
    
    # Ref seq and backtranseq backtranslation
    if on == "all":
        seq_df["Seq_similarity"] = seq_df.apply(lambda x: get_seq_similarity(x.Seq, translate_dna(x.Seq_dna)), axis=1) 
        avg_sim = np.mean(seq_df["Seq_similarity"].dropna())
        seq_ref = seq_df[seq_df["Reference_info"] == 1]
        seq_ref = seq_ref.drop_duplicates('Uniprot_ID')
        seq_backtr = seq_df[seq_df["Reference_info"] == 0]
        seq_backtr = seq_backtr.drop_duplicates('Uniprot_ID')
        avg_sim_ref = np.mean(seq_ref["Seq_similarity"].dropna())
        avg_sim_backtr = np.mean(seq_backtr["Seq_similarity"].dropna())
        if avg_sim < 1:                       
            logger.warning(f"Mismatch between protein and translated DNA. Overall average similatity: {avg_sim:.2f}")
        logger.debug(f"Average similarity ref coord backtranslation: {avg_sim_ref:.2f} for {len(seq_ref)} sequences.")      
        logger.debug(f"Average similarity backtranseq backtranslation: {avg_sim_backtr:.2f} for {len(seq_backtr)} sequences.") 
        
        return seq_df
         
    # Ref seq backtranslation only
    else:
        seq_df["Seq_similarity"] = seq_df.apply(lambda x: get_seq_similarity(x.Seq, translate_dna(x.Seq_dna)), axis=1)
        avg_sim = np.mean(seq_df["Seq_similarity"].dropna())
        if avg_sim < 1:
            tot_mismatch = sum(seq_df['Seq_similarity'].dropna() < 1)
            tot_seq = len(seq_df['Seq_similarity'].dropna())
            mismatch_ratio = tot_mismatch / tot_seq
            logger.warning(f"Mismatch between proteins and translated DNA by reference coordinates on {tot_mismatch} of {tot_seq}.")
            logger.warning(f"Protein and translated DNA mean similarity < 1: {avg_sim:.2f}")
            if mismatch_ratio > 0.4:
                logger.warning(f"Mismatch ratio > 0.4: {mismatch_ratio:.2f}. Dropping all {tot_mismatch} ref DNA sequnces...") 
                seq_df["Reference_info"] = 0
            else:
                logger.warning(f"Dropping {tot_mismatch} mismatching ref DNA sequnces...") 
                seq_df.loc[seq_df['Seq_similarity'] < 1, "Reference_info"] = 0
                
        return seq_df
                
                
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
    
    # Split by " "
    lst_extra_genes_rows = []
    for _, seq_row in seq_df.iterrows():
        uni_id = seq_row["Uniprot_ID"]
        gene_id = uniprot_to_gene_dict[uni_id]
        if pd.isnull(gene_id) == False:

            gene_id = gene_id.split(" ")
            if len(gene_id) > 1:
                for gene in gene_id:
                    if gene != seq_row["Gene"] and gene not in lst_added_genes:
                        
                        row = seq_row.copy()
                        row["Gene"] = gene
                        lst_extra_genes_rows.append(row)
                        lst_added_genes.append(gene)
    
    seq_df_extra_genes = pd.concat(lst_extra_genes_rows, axis=1).T
    
    # Split by "/"
    lst_extra_genes_rows = []
    for _, seq_row in seq_df.iterrows():
        uni_id = seq_row["Uniprot_ID"]
        gene_id = uniprot_to_gene_dict[uni_id]
        if pd.isnull(gene_id) == False:

            gene_id = gene_id.split("/")
            if len(gene_id) > 1:
                for gene in gene_id:
                    if gene != seq_row["Gene"] and gene not in lst_added_genes:
                        
                        row = seq_row.copy()
                        row["Gene"] = gene
                        lst_extra_genes_rows.append(row)
                        lst_added_genes.append(gene)
    
    seq_df_extra_genes2 = pd.concat(lst_extra_genes_rows, axis=1).T
    
    # Remove rows with multiple symbols and drop duplicated ones
    seq_df = pd.concat((seq_df, seq_df_extra_genes, seq_df_extra_genes2))
    seq_df = seq_df.dropna(subset=["Gene"])
    seq_df = seq_df[seq_df.apply(lambda x: len(x["Gene"].split(" ")), axis =1) == 1].reset_index(drop=True)
    seq_df = seq_df[seq_df.apply(lambda x: len(x["Gene"].split("/")), axis =1) == 1].reset_index(drop=True)
    seq_df = seq_df.drop_duplicates().reset_index(drop=True)
    
    return seq_df


def get_seq_df(input_dir, 
               output_seq_df, 
               uniprot_to_gene_dict = None, 
               organism = "Homo sapiens"):
    """
    Generate a dataframe including IDs mapping information (Gene and Uniprot_ID),
    AlphaFold 2 fragment number (F); the protein (Seq) and DNA sequence (DNA_seq) 
    obtained by EMBL backtranseq, as well as the genomic coordinate of the exons 
    (Chr and Exons_coord), which are used to compute the per-residue probability 
    of missense mutation. Also, for a subset of genes, it include the DNA sequence 
    of the reference genome (Seq_dna_ref) and its per-site trinucleotide context
    taing into account flanking regions at splicing sites (Tri_context).
    """

    # Load Uniprot ID to HUGO mapping
    uniprot_ids = os.listdir(input_dir)
    uniprot_ids = [uni_id.split("-")[1] for uni_id in list(set(uniprot_ids)) if ".pdb" in uni_id]
    if uniprot_to_gene_dict is not None and uniprot_to_gene_dict != "None":
        uniprot_to_gene_dict = json.load(open(uniprot_to_gene_dict)) 
    else:
        logger.debug("Retrieving Uniprot_ID to Hugo symbol mapping information..")
        uniprot_to_gene_dict = uniprot_to_hugo(uniprot_ids)  

    # Create a dataframe with protein sequences
    logger.debug("Generating sequence df...")
    seq_df = initialize_seq_df(input_dir, uniprot_to_gene_dict)
    
    # Add coordinates for mutability integration and ref DNA sequence extraction
    logger.debug("Retrieving exons coordinate..")
    coord_df = get_exons_coord(seq_df["Uniprot_ID"].unique())
    seq_df = seq_df.merge(coord_df, on=["Seq", "Uniprot_ID"], how="left").reset_index(drop=True)
    
    # Add ref DNA sequence and its per-site trinucleotide context
    logger.debug("Retrieving exons DNA sequence from reference genome..")
    if organism == "Homo sapiens":
        genome_fun = hg38
    elif organism == "Mus musculus":
        genome_fun = mm10
    else:
        raise RuntimeError(f"Failed to recognize '{organism}' as organism. Currently accepted ones are 'Homo sapiens' and 'Mus musculus'. Exiting...")
    seq_df = add_ref_dna_and_context(seq_df, genome_fun)
    
    # Check if ref DNA sequences match protein sequences
    seq_df = asssess_similarity(seq_df, on="ref")
    
    # Add DNA sequences for genes without available ref DNA using backtranseq
    logger.debug("Performing backtranseq backtranslation for genes without available ref DNA...")
    seq_df_backtranseq = seq_df[seq_df["Reference_info"] == 0].reset_index(drop=True)    
    seq_df_backtranseq = batch_backtranseq(seq_df_backtranseq, 500, organism=organism)    # TO DO: Check that the seq is overwritten
    seq_df_backtranseq["Tri_context"] = seq_df_backtranseq["Seq_dna"].apply(
        lambda x: ",".join(per_site_trinucleotide_context(x, no_flanks=True)))
    seq_df = pd.concat((seq_df[seq_df["Reference_info"] == 1], seq_df_backtranseq)
                       ).sort_values("Uniprot_ID").reset_index(drop=True)
    
    # Add multiple genes mapping to the same Uniprot_ID
    seq_df = add_extra_genes_to_seq_df(seq_df, uniprot_to_gene_dict)
    
    # Assess overal similarity and save
    seq_df = asssess_similarity(seq_df, on="all")
    seq_df.to_csv(output_seq_df, index=False)                    
    logger.debug(f"Sequences dataframe saved in: {output_seq_df}")
    