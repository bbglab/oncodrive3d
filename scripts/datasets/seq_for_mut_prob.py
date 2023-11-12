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

import requests
import sys
import numpy as np
import daiquiri
import pandas as pd
import requests
from progressbar import progressbar

from scripts import __logger_name__
from scripts.datasets.utils import (get_af_id_from_pdb,
                                    get_pdb_path_list_from_dir,
                                    get_seq_from_pdb, get_seq_similarity,
                                    translate_dna, uniprot_to_hugo)

logger = daiquiri.getLogger(__logger_name__ + ".build.seq_for_mut_prob")


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

    for path_structure in progressbar(list_prot_path):
        uniprot_id, f = get_af_id_from_pdb(path_structure).split("-F")
        gene = uniprot_to_gene_dict[uniprot_id]
        seq = "".join(list(get_seq_from_pdb(path_structure)))
        gene_lst.append(gene)
        uni_id_lst.append(uniprot_id)
        f_lst.append(f)
        seq_lst.append(seq)

    seq_df = pd.DataFrame({"Gene" : gene_lst, "Uniprot_ID" : uni_id_lst, "F" : f_lst, "Seq" : seq_lst}).sort_values(["Gene", "F"])

    return seq_df


def _uniprot_request_coord(lst_uniprot_ids):
    """
    Use Coordinates from EMBL-EBI Proteins API to get 
    a json including exons coordinate and protein info.
    
    https://www.ebi.ac.uk/proteins/api/doc/#coordinatesApi
    https://doi.org/10.1093/nar/gkx237
    """
    
    prot_request = [f"{prot}" if i == 0 else f"%2C{prot}" for i, prot in enumerate(lst_uniprot_ids)]
    requestURL = f"https://www.ebi.ac.uk/proteins/api/coordinates?offset=0&size=100&accession={''.join(prot_request)}"
    r = requests.get(requestURL, headers={ "Accept" : "application/json"})

    if not r.ok:
      r.raise_for_status()
      sys.exit()
    
    for dictio in json.loads(r.text):
        
        yield dictio


def get_batch_exons_coord(batch_ids):
    """
    Parse the json obtained from the Coordinates 
    service extracting exons coordinates and protein info.
    
    https://www.ebi.ac.uk/proteins/api/doc/#coordinatesApi
    https://doi.org/10.1093/nar/gkx237
    """
    
    lsts_gene = []
    lst_uni_id = []
    lst_seq = []
    lst_chr = []
    lst_ranges = []
    
    for dic in _uniprot_request_coord(batch_ids):
        gene = dic["gene"][0]["value"]
        uni_id = dic["accession"]
        seq = dic["sequence"]
        dic = dic["gnCoordinate"][0]["genomicLocation"]
        exons = dic["exon"]
        chromosome = dic["chromosome"]
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
            
        lsts_gene.append(gene)
        lst_uni_id.append(uni_id)
        lst_seq.append(seq)
        lst_chr.append(chromosome)
        lst_ranges.append(str(ranges))

    return pd.DataFrame({"Gene" : lsts_gene, "Uniprot_ID" : lst_uni_id, 
                         "Seq" : lst_seq, "Chr" : lst_chr, "Exons_coord" : lst_ranges})


def get_all_ids_coord(ids, batch_size=100):
    """
    Use the Coordinates service from Proteins API of EMBL-EBI to get 
    exons coordinate and proteins info of all provided Uniprot IDs.
    
    https://www.ebi.ac.uk/proteins/api/doc/#coordinatesApi
    https://doi.org/10.1093/nar/gkx237
    """
    
    lst_df = []
    batches_ids = [ids[i:i+batch_size] for i in range(0, len(ids), batch_size)]

    for batch_ids in progressbar(batches_ids):
        
        batch_df = get_batch_exons_coord(batch_ids)
        
        # Identify unmapped IDs and add them as NaN rows
        unmapped_ids = list(set(batch_ids).difference(set(batch_df.Uniprot_ID.unique())))
        
        nan_rows = pd.DataFrame({'Gene': np.repeat(np.nan, len(unmapped_ids)),
                                 'Uniprot_ID': unmapped_ids,
                                 'Seq': np.repeat(np.nan, len(unmapped_ids)),
                                 'Chr': np.repeat(np.nan, len(unmapped_ids)),
                                 'Exons_coord': np.repeat(np.nan, len(unmapped_ids))})

        batch_df = pd.concat([batch_df, nan_rows], ignore_index=True)
        lst_df.append(batch_df)
    
    return pd.concat(lst_df).reset_index(drop=True)


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
    response = requests.post(run_url, data=params)
    while str(response) != "<Response [200]>":
        time.sleep(5)
        response = requests.post(run_url, data=params)
    job_id = response.text.strip()

    # Wait for the job to complete
    status = ""
    while status != "FINISHED":
        time.sleep(5)
        result = requests.get(status_url + job_id)
        status = result.text.strip()

    # Retrieve the results of the job
    result = requests.get(result_url + job_id + "/out")
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
    for i, batch in progressbar(batches):
        
        # Avoid sending too many request
        if i+1 == 30:
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


def add_extra_genes_to_seq_df(seq_df, uniprot_to_gene_dict):
    """
    If multiple genes are mapping to a given Uniprot_ID, add 
    each gene name with corresponding sequence info to the seq_df.
    """
    
    lst_extra_genes = []

    for _, seq_row in seq_df.iterrows():
        uni_id = seq_row["Uniprot_ID"]
        gene_id = uniprot_to_gene_dict[uni_id]
        if pd.isnull(gene_id) == False:

            gene_id = gene_id.split(" ")
            if len(gene_id) > 1:
                for gene in gene_id:
                    if gene != seq_row["Gene"]:
                        lst_extra_genes.append((gene, seq_row["Uniprot_ID"], seq_row["F"], seq_row["Seq"], seq_row["Seq_dna"]))

    seq_df_extra_genes = pd.DataFrame(lst_extra_genes, columns=["Gene", "Uniprot_ID", "F", "Seq", "Seq_dna"])
    seq_df = pd.concat((seq_df, seq_df_extra_genes))
    
    # Remove rows with multiple symbles and drop duplicated ones
    seq_df = seq_df.dropna(subset=["Gene"])
    seq_df = seq_df[seq_df.apply(lambda x: len(x["Gene"].split(" ")), axis =1) == 1].reset_index(drop=True)
    seq_df = seq_df.drop_duplicates().reset_index(drop=True)
    
    return seq_df


def get_seq_df(input_dir, 
               output_seq_df, 
               uniprot_to_gene_dict = None, 
               organism = "Homo sapiens"):
    """
    Generate a dataframe including IDs mapping information (Gene and Uniprot_ID),
    AlphaFold 2 fragment number (F); the protein (Seq) and DNA sequence (DNA_seq) 
    as well as the genomic coordinate of the exons (Chr and Exons_coord), which 
    are used to compute the per-residue probability of missense mutation.
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
    
    # Annotate df with DNA sequences
    logger.debug("Performing back translation...")
    seq_df = batch_backtranseq(seq_df, 500, organism=organism)
    
    # Add multiple genes mapping to the same Uniprot_ID
    seq_df = add_extra_genes_to_seq_df(seq_df, uniprot_to_gene_dict)
    
    # Get coordinates for mutability integration (used to correct for seq depth in normal tissue)
    logger.debug("Retrieving exons coordinate..")
    coord_df = get_all_ids_coord(seq_df["Uniprot_ID"].unique())
    seq_df = seq_df.merge(coord_df.drop(columns=["Gene"]), on=["Seq", "Uniprot_ID"], how="left").reset_index(drop=True)
    
    # Save and assess similarity
    seq_df.to_csv(output_seq_df, index=False)
    sim_ratio = sum(seq_df.apply(lambda x: get_seq_similarity(x.Seq, translate_dna(x.Seq_dna)), axis=1)) / len(seq_df)
    if sim_ratio < 1:                       
        logger.warning(f"Some problems occurred during back translation: Protein and translated DNA similarity < 1: {sim_ratio}.")
    logger.debug(f"Dataframe including sequences is saved in: {output_seq_df}")