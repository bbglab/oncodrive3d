"""
Module to generate a pandas dataframe including identifiers 
mapped to protein and DNA sequences.

This dataframe is required to get the probability of each 
residue to mutate (missense mutation) based on the mutation 
profile (mutation rate in 96 trinucleotide contexts) of the 
cohort. The per-residue missense mutation probability of 
each protein is then used to get the probability of a 
certain volume to be hit by a missense mutation.
"""


import ast
import json
import logging
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import time

import daiquiri
import numpy as np
import pandas as pd
import requests
from bgreference import hg38, mm39
from Bio.Seq import Seq
from tqdm import tqdm

from scripts import __logger_name__
from scripts.datasets.utils import (
    download_single_file,
    get_af_id_from_pdb,
    get_pdb_path_list_from_dir,
    get_seq_from_pdb,
    uniprot_to_hugo,
)

logger = daiquiri.getLogger(__logger_name__ + ".build.seq_for_mut_prob")


_ENSEMBL_REST_SERVER = "https://rest.ensembl.org"
_ENSEMBL_REST_TIMEOUT = (10, 160)  # (connect, read) seconds
_ENSEMBL_REST_HEADERS = {"Accept": "text/plain"}
_ENSEMBL_CDS_MAX_CORES = 8


#===========
# region Initialize
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
            logger.debug(f"Request failed {e}: Retrying..")

    # Retrieve the results of the job
    status = "INIT"
    while status != "FINISHED":
        try:
            result = requests.get(result_url + job_id + "/out", timeout=160)
            status = "FINISHED"
        except requests.exceptions.RequestException as e:
            status = "ERROR"
            logger.debug(f"Request failed {e}: Retrying..")
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
                logger.debug(f"Error occurred after successfully sending request {r.raise_for_status()}: Retrying..")
                status = "ERROR"
        except requests.exceptions.RequestException as e:
            status = "ERROR"
            logger.debug(f"Request failed ({e}): Retrying..")

    for dictio in json.loads(r.text):

        yield dictio


def get_sorted_transcript_lst(dictio, canonical_ids):
    """
    Sort a list of tuple (index, transcript ID) placing the
    Ensembl canonical transcript as first element, if present.
    """

    lst_tuple_ix_id = [(i, coord_dict["ensemblTranscriptId"]) for (i, coord_dict) in enumerate(dictio)]

    return sorted(lst_tuple_ix_id, key=lambda x: x[1] in canonical_ids, reverse=True)


def get_batch_exons_coord(batch_ids, ens_canonical_transcripts_lst):
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

        # Iterate throug the transcripts (starting from Ensembl canonical one if present)
        sorted_transcript_lst = get_sorted_transcript_lst(dic["gnCoordinate"], ens_canonical_transcripts_lst)
        for i, ens_transcr_id in sorted_transcript_lst:

            ens_gene_id = dic["gnCoordinate"][i]["ensemblGeneId"]
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
            # If there is no transcript with matching sequence length, return NA
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


def get_exons_coord(ids, ens_canonical_transcripts_lst, batch_size=100):
    """
    Use the Coordinates service from Proteins API of EMBL-EBI to get
    exons coordinate and proteins info of all provided Uniprot IDs.

    https://www.ebi.ac.uk/proteins/api/doc/#coordinatesApi
    https://doi.org/10.1093/nar/gkx237
    """

    lst_df = []
    batches_ids = [ids[i:i+batch_size] for i in range(0, len(ids), batch_size)]

    for batch_ids in tqdm(batches_ids, total=len(batches_ids), desc="Adding exons coordinate"):

        batch_df = get_batch_exons_coord(batch_ids, ens_canonical_transcripts_lst)

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
            logger.warning(f"Error occurred during retrieving DNA seq from ref coordinates {transcript_id} {e}: Skipping..")

        return np.nan, np.nan, np.nan, -1


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
    seq_df = seq_df.merge(ref_dna_and_context, on=["Ens_Transcr_ID"], how="left").drop_duplicates().reset_index(drop=True)
    seq_df["Reference_info"] = seq_df["Reference_info"].fillna(-1).astype(int)
    seq_df.loc[seq_df["Reference_info"] == -1, "Ens_Transcr_ID"] = np.nan
    seq_df.loc[seq_df["Reference_info"] == -1, "Seq_dna"] = np.nan
    seq_df.loc[seq_df["Reference_info"] == -1, "Exons_coord"] = np.nan
    seq_df.loc[seq_df["Reference_info"] == -1, "Tri_context"] = np.nan

    return seq_df


#=========
# WRAPPERS
#=========

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

    if not lst_extra_genes_rows:
        return seq_df

    seq_df_extra_genes = pd.concat(lst_extra_genes_rows, axis=1).T

    # Remove rows with multiple symbols and drop duplicated ones
    seq_df = pd.concat((seq_df, seq_df_extra_genes))
    seq_df = seq_df.dropna(subset=["Gene"])
    seq_df = seq_df[seq_df.apply(lambda x: len(x["Gene"].split()), axis =1) == 1].reset_index(drop=True)
    seq_df = seq_df.drop_duplicates().reset_index(drop=True)

    return seq_df


def download_mane_summary(path_to_file, v=1.3, max_attempts=15, cores=1):
    """
    Download the summary.txt of the MANE release from NCBI.
    """

    mane_summary_url = f"https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_{v}/MANE.GRCh38.v{v}.summary.txt.gz"
    attempts = 0

    while not os.path.exists(path_to_file):
        download_single_file(mane_summary_url, path_to_file, threads=cores)
        attempts += 1
        if attempts >= max_attempts:
            raise RuntimeError(f"Failed to download MANE summary file after {max_attempts} attempts. Exiting..")
        time.sleep(5)


def select_uni_id(df, uniprot_ids):
    """
    For each row in `df`, normalize its 'Uniprot_ID' field—whether a single ID string
    or a semicolon-separated list—by selecting the first ID that exists in `available_ids`.
    If no IDs match, returns NaN for that row.
    """
    
    if sum(df.Uniprot_ID.str.split(";").apply(lambda x: len(x)) > 1) == 0:
        return df
        
    out = df.copy()
    
    def pick_first(uid_val):

        if pd.isna(uid_val):
            return np.nan
        
        # Split on ';' always: single IDs become a single-element list
        ids = str(uid_val).split(';')
        
        # Return first match or NaN
        return next((uid for uid in ids if uid in set(uniprot_ids)), np.nan)

    out['Uniprot_ID'] = out['Uniprot_ID'].apply(pick_first)
    return out


def load_custom_ens_prot_ids(path):
    """
    Read a CSV of custom PDB metadata and return a list of version‑stripped
    Ensembl Protein IDs.
    """
    if not os.path.isfile(path):
        logger.error(f"Custom PDB metadata path does not exist: {path!r}")
        raise FileNotFoundError(f"Custom PDB metadata not found: {path!r}")
    df = pd.read_csv(path)
    ids = (
        df["sequence"]
        .astype(str)
        .str.split(".", n=1)
        .str[0]
        .unique()
        .tolist()
        )
    
    return ids


def load_custom_symbol_map(path):
    """
    Read a samplesheet and return a mapping from ENSP (sequence) to gene symbol.
    Falls back to 'gene' column if 'symbol' is not present.
    """
    if not os.path.isfile(path):
        logger.error(f"Custom MANE metadata path does not exist: {path!r}")
        raise FileNotFoundError(f"Custom MANE metadata not found: {path!r}")
    df = pd.read_csv(path)
    if "sequence" not in df.columns:
        logger.debug("Custom MANE metadata missing 'sequence' column; skipping symbol mapping.")
        return {}

    symbol_col = None
    if "symbol" in df.columns:
        symbol_col = "symbol"
    elif "gene" in df.columns:
        symbol_col = "gene"
    else:
        logger.debug("Custom MANE metadata missing 'symbol'/'gene' column; skipping symbol mapping.")
        return {}

    seq = df["sequence"].astype(str).str.split(".", n=1).str[0]
    sym = df[symbol_col]
    mapping = {}
    for k, v in zip(seq, sym):
        if pd.isna(v):
            continue
        v = str(v).strip()
        if not v:
            continue
        mapping[k] = v
    return mapping


def get_mane_to_af_mapping(
    datasets_dir, 
    uniprot_ids, 
    include_not_af=False, 
    mane_version=1.4,
    custom_mane_metadata_path=None, 
    cores=1
    ):
    """
    Get a dataframe to map genes, MANE transcript IDs, and AlphaFold structures.
    If custom_mane_metadata_path is provided, the Ensembl protein IDs within the 
    dataframe will be used to overwrite the `Uniprot_ID` column (version stripped).

    Parameters
    ----------
    datasets_dir : str
        Root folder where MANE and UniProt files live.
    uniprot_ids : set
        UniProt accessions already downloaded locally (used to pick among multiples).
    include_not_af : bool, default False
        If True, also return MANE entries without AlphaFold models.
    mane_version : float, default 1.4
        Which MANE release to fetch.
    custom_mane_metadata_path : str, optional
        If provided, these Ensembl protein IDs will be used to overwrite
        the `Uniprot_ID` column (version stripped) for those entries.

    Returns
    -------
    pandas.DataFrame or (DataFrame, DataFrame)
        The merged mapping (and, if `include_not_af`, a second DF of missing entries in AF).
    """

    mane_to_af = pd.read_csv(os.path.join(datasets_dir, "mane_refseq_prot_to_alphafold.csv"))
    mane_to_af = mane_to_af.rename(columns={"refseq_prot" : "Refseq_prot",
                                            "uniprot_accession" : "Uniprot_ID"}).drop(columns=["alphafold"])
    path_mane_summary = os.path.join(datasets_dir, "mane_summary.txt.gz")
    if not os.path.exists(path_mane_summary):
        download_mane_summary(path_mane_summary, mane_version, cores=cores)

    mane_summary = pd.read_csv(path_mane_summary, compression='gzip', sep="\t")
    
    mane_summary = mane_summary.rename(columns={
        "symbol" : "Gene",
        "RefSeq_prot" : "Refseq_prot",
        "Ensembl_Gene" : "Ens_Gene_ID",
        "Ensembl_prot" : "Ens_Prot_ID",
        "Ensembl_nuc" : "Ens_Transcr_ID",
        "GRCh38_chr" : "Chr",
        "chr_strand" : "Reverse_strand"
        })
    mane_summary = mane_summary[[
        "Gene", 
        "HGNC_ID", 
        "Ens_Gene_ID", 
        "Ens_Prot_ID", 
        "Refseq_prot", 
        "Ens_Transcr_ID", 
        "Chr", 
        "Reverse_strand"
        ]]
    
    mane_mapping = mane_summary.merge(mane_to_af, how="left", on="Refseq_prot")
    mane_mapping.Reverse_strand = mane_mapping.Reverse_strand.map({"+" : 0, "-" : 1})
    mane_mapping.Ens_Gene_ID = mane_mapping.Ens_Gene_ID.apply(lambda x: x.split(".")[0])
    mane_mapping.Ens_Transcr_ID = mane_mapping.Ens_Transcr_ID.apply(lambda x: x.split(".")[0])

    # Override Uniprot_ID with Ens_Prot_ID for the custom PDB structures
    if custom_mane_metadata_path is not None:
        custom_ids = load_custom_ens_prot_ids(custom_mane_metadata_path)
        base_ens = mane_mapping["Ens_Prot_ID"].str.split(".", n=1).str[0]
        mask = base_ens.isin(custom_ids)
        mane_mapping.loc[mask, "Uniprot_ID"] = base_ens[mask]
        if logger.isEnabledFor(logging.DEBUG):
            summary_ens = set(mane_summary["Ens_Prot_ID"].astype(str).str.split(".", n=1).str[0])
            missing_custom = sorted(set(custom_ids) - summary_ens)
            if missing_custom:
                preview = ", ".join(missing_custom[:10])
                suffix = "..." if len(missing_custom) > 10 else ""
                logger.debug(
                    "Custom MANE ENSP IDs not found in MANE summary (%s): %s%s",
                    len(missing_custom),
                    preview,
                    suffix,
                )

    # Select available Uniprot ID, fist one if multiple are present
    mane_mapping = mane_mapping.dropna(subset=["Uniprot_ID"]).reset_index(drop=True)
    mane_mapping = select_uni_id(mane_mapping, uniprot_ids)
    mane_mapping = mane_mapping[mane_mapping["Uniprot_ID"].isin(uniprot_ids)]

    # Also return a dataframe with entries not in AF
    if include_not_af:
        mane_not_af = mane_summary[~mane_summary.Gene.isin(mane_mapping.Gene)].reset_index(drop=True)
        return mane_mapping, mane_not_af

    else:
        return mane_mapping


def download_biomart_metadata(path_to_file, max_attempts=5, wait_seconds=10):
    """
    Query biomart to get the list of transcript corresponding to the downloaded
    structures (a few structures are missing) and other information.
    """

    base_archive = "http://jan2024.archive.ensembl.org"
    base_latest = "https://www.ensembl.org"
    query = (
        '/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?>'
        '<!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="0" count="" '
        'datasetConfigVersion="0.6"><Dataset name="hsapiens_gene_ensembl" interface="default">'
        '<Attribute name="ensembl_gene_id"/><Attribute name="ensembl_transcript_id"/>'
        '<Attribute name="transcript_is_canonical"/><Attribute name="external_gene_name"/>'
        '<Attribute name="external_gene_source"/><Attribute name="hgnc_id"/>'
        '<Attribute name="uniprot_gn_id"/><Attribute name="uniprotswissprot"/>'
        '<Attribute name="external_synonym"/></Dataset></Query>'
    )
    url = f"{base_archive}{query}"
    fallback_url = f"{base_latest}{query}"
    logger.debug("Starting BioMart metadata download to %s (archive: %s, latest: %s).", path_to_file, base_archive, base_latest)

    if shutil.which("wget") is None:
        logger.warning("wget not found; falling back to Python downloader for BioMart metadata.")
        last_exc = None
        ssl_verify_archive = url.startswith("https://")
        for attempt in range(1, max_attempts + 1):
            logger.debug("Starting BioMart download attempt %s/%s (archive).", attempt, max_attempts)
            try:
                download_single_file(url, path_to_file, threads=4, ssl=ssl_verify_archive)
                return
            except Exception as exc:
                last_exc = exc
                logger.warning(
                    "BioMart download failed (attempt %s/%s). Retrying in %ss... Error: %s",
                    attempt,
                    max_attempts,
                    wait_seconds,
                    exc,
                )
                logger.debug("BioMart download exception details:", exc_info=True)
                time.sleep(wait_seconds)

        logger.warning("Falling back to latest Ensembl BioMart URL after failure on %s.", base_archive)
        if os.path.exists(path_to_file):
            try:
                os.remove(path_to_file)
            except OSError as exc:
                logger.warning(
                    "Failed to remove partial BioMart metadata file %s before fallback: %s",
                    path_to_file,
                    exc,
                )
        ssl_verify_fallback = fallback_url.startswith("https://")
        for attempt in range(1, max_attempts + 1):
            logger.debug("Starting BioMart download attempt %s/%s (latest).", attempt, max_attempts)
            try:
                download_single_file(fallback_url, path_to_file, threads=4, ssl=ssl_verify_fallback)
                return
            except Exception as exc:
                last_exc = exc
                logger.warning(
                    "Fallback BioMart download failed (attempt %s/%s). Retrying in %ss... Error: %s",
                    attempt,
                    max_attempts,
                    wait_seconds,
                    exc,
                )
                logger.debug("Fallback BioMart download exception details:", exc_info=True)
                time.sleep(wait_seconds)

        raise RuntimeError(
            f"Failed to download BioMart metadata after {max_attempts} attempts on archive and "
            f"{max_attempts} attempts on latest."
        ) from last_exc

    command = [
        "wget",
        "--no-hsts",
        "--continue",
        "--read-timeout=120",
        "--timeout=120",
        "--tries=1",
        "-O",
        path_to_file,
        url,
    ]

    for attempt in range(1, max_attempts + 1):
        logger.debug("Starting BioMart wget attempt %s/%s (archive).", attempt, max_attempts)
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode == 0:
            return
        stderr = (result.stderr or "").strip()
        if stderr:
            logger.warning(
                "BioMart download failed (attempt %s/%s, return code %s). stderr: %s",
                attempt,
                max_attempts,
                result.returncode,
                stderr,
            )
        else:
            logger.warning(
                "BioMart download failed (attempt %s/%s, return code %s). Retrying in %ss...",
                attempt,
                max_attempts,
                result.returncode,
                wait_seconds,
            )
        if result.stdout:
            logger.debug("BioMart wget stdout (attempt %s/%s): %s", attempt, max_attempts, result.stdout.strip())
        time.sleep(wait_seconds)

    logger.warning("Falling back to latest Ensembl BioMart URL after failure on %s.", base_archive)
    if os.path.exists(path_to_file):
        try:
            os.remove(path_to_file)
        except OSError as exc:
            logger.warning(
                "Failed to remove partial BioMart metadata file %s before fallback: %s",
                path_to_file,
                exc,
            )
    command[-1] = fallback_url
    for attempt in range(1, max_attempts + 1):
        logger.debug("Starting BioMart wget attempt %s/%s (latest).", attempt, max_attempts)
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode == 0:
            return
        stderr = (result.stderr or "").strip()
        if stderr:
            logger.warning(
                "Fallback BioMart download failed (attempt %s/%s, return code %s). stderr: %s",
                attempt,
                max_attempts,
                result.returncode,
                stderr,
            )
        else:
            logger.warning(
                "Fallback BioMart download failed (attempt %s/%s, return code %s). Retrying in %ss...",
                attempt,
                max_attempts,
                result.returncode,
                wait_seconds,
            )
        if result.stdout:
            logger.debug(
                "Fallback BioMart wget stdout (attempt %s/%s): %s",
                attempt,
                max_attempts,
                result.stdout.strip(),
            )
        time.sleep(wait_seconds)

    raise RuntimeError(
        f"Failed to download BioMart metadata after {max_attempts} attempts on archive and "
        f"{max_attempts} attempts on latest."
    )


def get_biomart_metadata(datasets_dir, uniprot_ids):
    """
    Download a dataframe including ensembl canonical transcript IDs,
    HGNC IDs, Uniprot IDs, and other useful information.
    """

    try:
        path_biomart_metadata = os.path.join(datasets_dir, "biomart_metadata.tsv")
        if not os.path.exists(path_biomart_metadata):
            download_biomart_metadata(path_biomart_metadata)

        # Parse
        biomart_df = pd.read_csv(path_biomart_metadata, sep="\t", header=None, low_memory=False)
        biomart_df.columns = ["Ens_Gene_ID", "Ens_Transcr_ID", "Ens_Canonical", "Gene", "Gene_source", "HGNC_ID", "Uniprot_ID", "UniprotKB_ID", "Gene_synonym"]
        biomart_df = biomart_df.dropna(subset=["Uniprot_ID", "UniprotKB_ID"], how='all')
        biomart_df = biomart_df[biomart_df["Gene_source"] == "HGNC Symbol"].drop(columns=["Gene_source"])
        biomart_df.reset_index(inplace=True, drop=True)

        # Filter
        uniprot_ids = pd.Series(uniprot_ids)
        uniprot_kb_ids = uniprot_ids[(uniprot_ids.isin(biomart_df.UniprotKB_ID))]
        uniprot_notkb_ids = uniprot_ids[~(uniprot_ids.isin(biomart_df.UniprotKB_ID)) & (uniprot_ids.isin(biomart_df.Uniprot_ID))]
        biomart_uniprotkb = biomart_df[biomart_df.UniprotKB_ID.isin(uniprot_kb_ids)].drop(
            columns=["Uniprot_ID"]).rename(columns={"UniprotKB_ID" : "Uniprot_ID"}).drop_duplicates()
        biomart_uniprot = biomart_df[biomart_df.Uniprot_ID.isin(uniprot_notkb_ids)].drop(columns=["UniprotKB_ID"]).drop_duplicates()
        biomart_df = pd.concat((biomart_uniprotkb, biomart_uniprot)).reset_index(drop=True)

        # Output
        biomart_df.to_csv(path_biomart_metadata, sep="\t", index=False)
        canonical_transcripts = biomart_df.Ens_Transcr_ID[biomart_df["Ens_Canonical"] == 1].unique()

    except Exception as e:
        logger.warning(
            "Metadata from BioMart could not be downloaded; transcript IDs will not be prioritized by canonical.  Error was: %s",
            e,
            exc_info=True
        )
        canonical_transcripts = []

    return canonical_transcripts


def get_ref_dna_from_ensembl_batch(transcript_ids, max_attempts=8, wait_seconds=3):
    """
    Retrieve Ensembl CDS DNA sequences for up to 50 stable IDs in a single request.

    Ensembl REST docs: POST /sequence/id (max POST size = 50).
    """

    pid = os.getpid()
    start_time = time.perf_counter()

    if not transcript_ids:
        return []

    # Keep output aligned with the input order (including any NA values).
    results = [np.nan] * len(transcript_ids)
    request_ids = []
    request_positions = []
    for pos, tid in enumerate(transcript_ids):
        if pd.isna(tid):
            continue
        tid_str = str(tid)
        request_ids.append(tid_str)
        request_positions.append(pos)

    if not request_ids:
        return results

    if len(request_ids) > 50:
        raise ValueError(f"Ensembl POST /sequence/id supports max 50 IDs per request; got {len(request_ids)}.")

    url = f"{_ENSEMBL_REST_SERVER}/sequence/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    payload = {"ids": request_ids}
    params = {"type": "cds"}

    last_error = None
    for attempt in range(1, max_attempts + 1):
        try:
            r = requests.post(url, headers=headers, json=payload, params=params, timeout=_ENSEMBL_REST_TIMEOUT)
            if r.status_code == 429:
                retry_after = r.headers.get("Retry-After")
                try:
                    retry_after = float(retry_after)
                except (TypeError, ValueError):
                    retry_after = wait_seconds
                if attempt >= max_attempts:
                    logger.warning(
                        "Ensembl CDS batch rate limited (pid=%s, attempt=%s/%s). Giving up after %ss.",
                        pid,
                        attempt,
                        max_attempts,
                        retry_after,
                    )
                    return results
                logger.warning(
                    "Ensembl CDS batch rate limited (pid=%s, attempt=%s/%s). Retrying after %ss.",
                    pid,
                    attempt,
                    max_attempts,
                    retry_after,
                )
                time.sleep(retry_after)
                continue

            r.raise_for_status()
            decoded = r.json()

            seq_by_id = {}
            errors_by_id = {}
            if isinstance(decoded, dict):
                decoded = [decoded]
            if not isinstance(decoded, list):
                raise ValueError(f"Unexpected Ensembl response type: {type(decoded).__name__}")

            for item in decoded:
                if not isinstance(item, dict):
                    continue

                item_id = item.get("id")
                item_query = item.get("query")
                item_keys = [k for k in (item_query, item_id) if k]
                if not item_keys:
                    continue

                if "seq" in item and item.get("seq") is not None:
                    seq_val = item.get("seq")
                    for key in item_keys:
                        seq_by_id[key] = seq_val
                elif "error" in item:
                    err_val = item.get("error")
                    for key in item_keys:
                        errors_by_id[key] = err_val

            missing = 0
            for pos, tid in zip(request_positions, request_ids):
                seq_dna = seq_by_id.get(tid)
                if seq_dna is None and "." in tid:
                    seq_dna = seq_by_id.get(tid.split(".", 1)[0])
                if not seq_dna:
                    missing += 1
                    continue
                results[pos] = seq_dna[:-3] if len(seq_dna) >= 3 else np.nan

            elapsed = time.perf_counter() - start_time
            if missing > 0 and logger.isEnabledFor(logging.DEBUG):
                example_missing = [tid for tid in request_ids if tid not in seq_by_id][:5]
                logger.debug(
                    "Ensembl CDS batch completed with missing IDs (pid=%s, elapsed=%.2fs, requested=%s, missing=%s). Example missing: %s",
                    pid,
                    elapsed,
                    len(request_ids),
                    missing,
                    example_missing,
                )
                if errors_by_id:
                    example_errors = {k: errors_by_id[k] for k in list(errors_by_id)[:3]}
                    logger.debug("Ensembl CDS batch errors (pid=%s): %s", pid, example_errors)
            else:
                logger.debug(
                    "Ensembl CDS batch completed (pid=%s, elapsed=%.2fs, requested=%s)",
                    pid,
                    elapsed,
                    len(request_ids),
                )

            return results

        except (requests.exceptions.RequestException, ValueError, json.JSONDecodeError) as exc:
            last_error = exc
            if attempt < max_attempts:
                logger.debug(
                    "Ensembl CDS batch request failed (pid=%s, attempt=%s/%s): %s. Retrying in %ss...",
                    pid,
                    attempt,
                    max_attempts,
                    exc,
                    wait_seconds,
                )
                time.sleep(wait_seconds)
                continue

            elapsed = time.perf_counter() - start_time
            logger.warning(
                "Ensembl CDS batch failed (pid=%s, elapsed=%.2fs, requested=%s). Last error: %s",
                pid,
                elapsed,
                len(request_ids),
                last_error,
            )

    return results


def get_ref_dna_from_ensembl(transcript_id):
    """
    Use Ensembl GET sequence rest API to obtain CDS DNA
    sequence from Ensembl transcript ID.

    https://rest.ensembl.org/documentation/info/sequence_id
    """

    pid = os.getpid()
    start_time = time.perf_counter()
    failures = 0

    if pd.isna(transcript_id):
        logger.debug("Ensembl CDS start: <NA> (pid=%s) -> skipping", pid)
        return np.nan

    transcript_id = str(transcript_id)
    logger.debug("Ensembl CDS start: %s (pid=%s)", transcript_id, pid)

    url = f"{_ENSEMBL_REST_SERVER}/sequence/id/{transcript_id}?type=cds"

    status = "INIT"
    last_error = None
    while status != "FINISHED":

        try:
            r = requests.get(url, headers=_ENSEMBL_REST_HEADERS, timeout=_ENSEMBL_REST_TIMEOUT)
            if not r.ok:
                r.raise_for_status()

                status = "ERROR"
                sys.exit()
            else:
                status = "FINISHED"

        except requests.exceptions.RequestException as e:
            failures += 1
            last_error = e
            status = "ERROR"
            if failures % 10 == 0:
                logger.debug(
                    "Failed to retrieve sequence for %s (pid=%s, failures=%s) %s: Retrying..",
                    transcript_id,
                    pid,
                    failures,
                    e,
                )
                time.sleep(5)
            if failures == 100:
                elapsed = time.perf_counter() - start_time
                logger.warning(
                    "Ensembl CDS failed: %s (pid=%s, elapsed=%.2fs, failures=%s). Last error: %s",
                    transcript_id,
                    pid,
                    elapsed,
                    failures,
                    last_error,
                )
                return np.nan

            time.sleep(1)

    text = r.text.strip()
    if text.startswith(">"):
        seq_dna = "".join(text.splitlines()[1:])
    else:
        seq_dna = "".join(text.splitlines())

    elapsed = time.perf_counter() - start_time
    if failures > 0:
        logger.info(
            "Ensembl CDS completed: %s (pid=%s, elapsed=%.2fs, failures=%s)",
            transcript_id,
            pid,
            elapsed,
            failures,
        )
    else:
        logger.debug(
            "Ensembl CDS completed: %s (pid=%s, elapsed=%.2fs, failures=%s)",
            transcript_id,
            pid,
            elapsed,
            failures,
        )

    return seq_dna[:-3] if len(seq_dna) >= 3 else np.nan


def get_ref_dna_from_ensembl_mp(seq_df, cores):
    """
    Multiple processing function to use Ensembl GET sequence
    rest API to obtain CDS DNA sequence from Ensembl transcript ID.

    https://rest.ensembl.org/documentation/info/sequence_id
    """

    seq_df = seq_df.copy()
    transcript_ids = seq_df.Ens_Transcr_ID.tolist()
    total = len(transcript_ids)
    if total == 0:
        seq_df["Seq_dna"] = []
        return seq_df

    if cores > _ENSEMBL_CDS_MAX_CORES:
        logger.info(
            "Capping Ensembl CDS batch workers from %s to %s.",
            cores,
            _ENSEMBL_CDS_MAX_CORES,
        )
        cores = _ENSEMBL_CDS_MAX_CORES

    logger.debug("Retrieving CDS DNA from Ensembl for %s transcripts with %s cores.", total, cores)

    if cores <= 1:
        results = []
        batch_size = 50
        with tqdm(total=total, desc="Ensembl CDS") as pbar:
            for i in range(0, total, batch_size):
                batch_ids = transcript_ids[i : i + batch_size]
                batch_results = get_ref_dna_from_ensembl_batch(batch_ids)
                results.extend(batch_results)
                pbar.update(len(batch_ids))
        seq_df["Seq_dna"] = results
        logger.debug("Completed Ensembl CDS retrieval.")
        return seq_df

    batch_size = 50
    batches = [transcript_ids[i : i + batch_size] for i in range(0, total, batch_size)]
    results = []
    with multiprocessing.Pool(processes=cores) as pool:
        results_iter = pool.imap(get_ref_dna_from_ensembl_batch, batches)
        with tqdm(total=total, desc="Ensembl CDS") as pbar:
            for batch_ids, batch_results in zip(batches, results_iter):
                results.extend(batch_results)
                pbar.update(len(batch_ids))
    seq_df["Seq_dna"] = results

    logger.debug("Completed Ensembl CDS retrieval.")
    return seq_df


def drop_gene_duplicates(df):
    """
    Issue: multiple Uniprot IDs might map to the same HUGO symbol.
    Drop Uniprot ID of gene duplicates by prioritizing reference
    info status (Proteins API > Ensembl API > Backtranseq API).
    """

    df = df.copy()
    df.Uniprot_ID = df.Uniprot_ID.str.replace(";", "")
    df.Gene = df.Gene.str.replace(";", "")
    df = df.sort_values(["Gene", "Reference_info"], ascending=False).drop_duplicates(subset='Gene').reset_index(drop=True)

    # Check if there are still duplicates
    n_duplicates = sum(df["Gene"].value_counts() > 1)
    if n_duplicates > 0:
        logger.warning(f"Found {n_duplicates} duplicates gene entries: Mapping HUGO Symbol to protein info might be affected.")
    else:
        logger.debug("Duplicates gene entries correctly removed!")

    return df


def mane_uniprot_to_hugo(uniprot_ids, mane_mapping):
    """
    Generate a mapping from Uniprot IDs to Hugo Symbols for MANE 
    transcripts associated structures.
    """
    
    sel = mane_mapping.loc[mane_mapping["Uniprot_ID"].isin(uniprot_ids),["Uniprot_ID", "Gene"]]
    return dict(zip(sel["Uniprot_ID"], sel["Gene"]))


def process_seq_df(seq_df,
                   datasets_dir,
                   organism,
                   uniprot_to_gene_dict):
    """
    Retrieve DNA sequence and tri-nucleotide context
    for each structure in the initialized dataframe
    prioritizing structures obtained from transcripts
    whose exon coordinates are available in the Proteins API.

    Canonical transcript metadata is retrieved internally from BioMart
    for the provided dataset directory and Uniprot IDs.

    Reference_info labels:
        1 : Transcript ID, exons coord, seq DNA obtained from Proteins API
       -1 : Not available transcripts, seq DNA retrieved from Backtranseq API
    """

    if seq_df.empty:
        logger.error("No sequences to process in process_seq_df; this should not happen.")
        raise RuntimeError("Empty sequence dataframe: no structures to process.")

    ens_canonical_transcripts_lst = get_biomart_metadata(datasets_dir, seq_df["Uniprot_ID"].unique())

    # Process entries in Proteins API (Reference_info 1)
    #---------------------------------------------------

    # Add coordinates for mutability integration (entries in Proteins API)
    logger.debug(f"Retrieving CDS DNA seq from reference genome (Proteins API): {len(seq_df['Uniprot_ID'].unique())} structures..")
    coord_df = get_exons_coord(seq_df["Uniprot_ID"].unique(), ens_canonical_transcripts_lst)
    seq_df = seq_df.merge(coord_df, on=["Seq", "Uniprot_ID"], how="left").reset_index(drop=True)

    # Add ref DNA seq and its per-site trinucleotide context (entries in Proteins API)
    if organism == "Homo sapiens":
        logger.debug("Loading reference genome hg38..")
        genome_fun = hg38
    elif organism == "Mus musculus":
        logger.debug("Loading reference genome mm39..")
        genome_fun = mm39
    else:
        raise RuntimeError(f"Failed to recognize '{organism}' as organism. Currently accepted ones are 'Homo sapiens' and 'Mus musculus'. Exiting..")
    seq_df = add_ref_dna_and_context(seq_df, genome_fun)
    seq_df_uniprot = seq_df[seq_df["Reference_info"] == 1]
    seq_df_not_uniprot = seq_df[seq_df["Reference_info"] == -1]


    # Process entries not in Proteins API (Reference_info -1)
    #------------------------------------------------------------

    if len(seq_df_not_uniprot) > 0:
        # Add DNA seq from Backtranseq for any other entry
        logger.debug(f"Retrieving CDS DNA seq for entries without available transcript ID (Backtranseq API): {len(seq_df_not_uniprot)} structures..")
        seq_df_not_uniprot = batch_backtranseq(seq_df_not_uniprot, 500, organism=organism)

        # Get trinucleotide context
        seq_df_not_uniprot["Tri_context"] = np.nan
        valid_seq_mask = seq_df_not_uniprot["Seq_dna"].notna()
        seq_df_not_uniprot.loc[valid_seq_mask, "Tri_context"] = seq_df_not_uniprot.loc[valid_seq_mask, "Seq_dna"].apply(
            lambda x: ",".join(per_site_trinucleotide_context(x, no_flanks=True)))


    # Prepare final output
    #---------------------

    # Concat the dfs, expand multiple genes associated to the same structure, keep only one structure for each gene
    seq_df = pd.concat((seq_df_uniprot, seq_df_not_uniprot)).reset_index(drop=True)
    logger_report = ", ".join([f"{v}: {c}" for (v, c) in zip(seq_df.Reference_info.value_counts().index,
                                                              seq_df.Reference_info.value_counts().values)])
    logger.info(f"Built of sequence dataframe completed. Retrieved {len(seq_df)} structures ({logger_report})")
    seq_df = add_extra_genes_to_seq_df(seq_df, uniprot_to_gene_dict)
    pre_drop = len(seq_df)
    seq_df = drop_gene_duplicates(seq_df)
    logger.debug(
        "Duplicate gene removal: %s removed (from %s to %s).",
        pre_drop - len(seq_df),
        pre_drop,
        len(seq_df),
    )

    return seq_df


def process_seq_df_mane(seq_df,
                        datasets_dir,
                        uniprot_to_gene_dict,
                        mane_mapping,
                        mane_mapping_not_af,
                        mane_only=False,
                        num_cores=1):
    """
    Retrieve DNA sequence and tri-nucleotide context
    for each structure in the initialized dataframe
    prioritizing MANE associated structures and metadata.

    Reference_info labels:
        1  : Transcript ID, exons coord, seq DNA obtained from Proteins API
        0  : Transcript ID retrieved from MANE and seq DNA from Ensembl
        -1 : Not available transcripts, seq DNA retrieved from Backtranseq API
    """

    seq_df_mane = seq_df[seq_df.Uniprot_ID.isin(mane_mapping.Uniprot_ID)].reset_index(drop=True)
    seq_df_nomane = seq_df[~seq_df.Uniprot_ID.isin(mane_mapping.Uniprot_ID)].reset_index(drop=True)

    # Seq df MANE
    # -----------
    seq_df_mane = seq_df_mane.drop(columns=["Gene"]).merge(mane_mapping, how="left", on="Uniprot_ID")
    seq_df_mane["Reference_info"] = 0

    if seq_df_mane.empty:
        logger.warning("No MANE sequences to process; skipping Ensembl CDS retrieval.")
        if "Seq_dna" not in seq_df_mane.columns:
            seq_df_mane["Seq_dna"] = pd.Series(dtype=object)
    else:
        # Add DNA seq from Ensembl for structures with available transcript ID
        logger.debug(f"Retrieving CDS DNA seq from transcript ID (Ensembl API): {len(seq_df_mane)} structures..")
        seq_df_mane = get_ref_dna_from_ensembl_mp(seq_df_mane, cores=num_cores)

        # Retry missing entries one-by-one using single-request API
        missing_mask = seq_df_mane["Seq_dna"].isna()
        if missing_mask.any():
            logger.debug(
                "Retrying %s missing Ensembl CDS entries one by one.",
                int(missing_mask.sum()),
            )
            recovered = 0
            for idx in seq_df_mane[missing_mask].index:
                seq_dna = get_ref_dna_from_ensembl(seq_df_mane.at[idx, "Ens_Transcr_ID"])
                if pd.notna(seq_dna):
                    recovered += 1
                seq_df_mane.at[idx, "Seq_dna"] = seq_dna
            if recovered > 0:
                logger.debug("Recovered %s missing Ensembl CDS entries after single-request retry.", recovered)

        # Set failed and len-mismatching entries as no-transcripts entries
        seq_len = seq_df_mane["Seq"].str.len()
        dna_len = seq_df_mane["Seq_dna"].str.len()
        failed_nan = seq_df_mane["Seq_dna"].isna()
        failed_mismatch = (~failed_nan) & (dna_len / 3 != seq_len)
        failed_ix = failed_nan | failed_mismatch
        logger.debug(
            "Ensembl CDS failures: total=%s (missing=%s, length_mismatch=%s).",
            int(failed_ix.sum()),
            int(failed_nan.sum()),
            int(failed_mismatch.sum()),
        )
        if sum(failed_ix) > 0:
            seq_df_mane_failed = seq_df_mane[failed_ix]
            seq_df_mane = seq_df_mane[~failed_ix]
            seq_df_mane_failed = seq_df_mane_failed.drop(columns=[
                "Ens_Gene_ID",  
                "Ens_Transcr_ID", 
                "Reverse_strand",
                "Chr", 
                "Refseq_prot", 
                "Reference_info", 
                "Seq_dna"
                ])
            seq_df_nomane = pd.concat((seq_df_nomane, seq_df_mane_failed))
            logger.debug(
                "Moved %s failed MANE entries to non-MANE pool.",
                len(seq_df_mane_failed),
            )

    # Seq df not MANE
    # ---------------
    if seq_df_nomane.empty:
        logger.debug("No non-MANE sequences to process; skipping Proteins/Backtranseq retrieval.")
        seq_df_nomane_tr = seq_df_nomane.copy()
        seq_df_nomane_notr = seq_df_nomane.copy()
    else:
        before_nomane = len(seq_df_nomane)
        seq_df_nomane = add_extra_genes_to_seq_df(seq_df_nomane, uniprot_to_gene_dict)         # Filter out genes with NA
        after_extra = len(seq_df_nomane)
        if not mane_only:
            logger.debug("Filtering non-MANE entries using mane_mapping_not_af (gene whitelist).")
            seq_df_nomane = seq_df_nomane[seq_df_nomane.Gene.isin(mane_mapping_not_af.Gene)]       # Filter out genes that are not in MANE list
            logger.debug(
                "Non-MANE pool sizes: initial=%s, after_extra_genes=%s, after_mane_filter=%s.",
                before_nomane,
                after_extra,
                len(seq_df_nomane),
            )
        else:
            logger.debug(
                "Non-MANE pool sizes: initial=%s, after_extra_genes=%s (mane_only, no MANE filter applied).",
                before_nomane,
                after_extra,
            )

        if seq_df_nomane.empty:
            logger.debug("No non-MANE sequences after filtering; skipping Proteins/Backtranseq retrieval.")
            seq_df_nomane_tr = seq_df_nomane.copy()
            seq_df_nomane_notr = seq_df_nomane.copy()
        else:
            ens_canonical_transcripts_lst = get_biomart_metadata(datasets_dir, seq_df_nomane["Uniprot_ID"].unique())
            # Retrieve seq from coordinates
            logger.debug(f"Retrieving CDS DNA seq from reference genome (Proteins API): {len(seq_df_nomane['Uniprot_ID'].unique())} structures..")
            coord_df = get_exons_coord(seq_df_nomane["Uniprot_ID"].unique(), ens_canonical_transcripts_lst)
            seq_df_nomane = seq_df_nomane.merge(coord_df, on=["Seq", "Uniprot_ID"], how="left").reset_index(drop=True)  # Discard entries whose Seq obtained by Proteins API don't exactly match the one in structure
            seq_df_nomane = add_ref_dna_and_context(seq_df_nomane, hg38)
            seq_df_nomane_tr = seq_df_nomane[seq_df_nomane["Reference_info"] == 1]
            seq_df_nomane_notr = seq_df_nomane[seq_df_nomane["Reference_info"] == -1]

            # Add DNA seq from Backtranseq for any other entry
            if len(seq_df_nomane_notr) > 0:
                logger.debug(f"Retrieving CDS DNA seq for genes without available transcript ID (Backtranseq API): {len(seq_df_nomane_notr)} structures..")
                seq_df_nomane_notr = batch_backtranseq(seq_df_nomane_notr, 500, organism="Homo sapiens")

    # Get trinucleotide context
    seq_df_not_uniprot = pd.concat((seq_df_mane, seq_df_nomane_notr))
    if "Seq_dna" not in seq_df_not_uniprot.columns:
        seq_df_not_uniprot["Seq_dna"] = pd.Series(dtype=object)
    seq_df_not_uniprot["Tri_context"] = np.nan
    valid_seq_mask = seq_df_not_uniprot["Seq_dna"].notna()
    seq_df_not_uniprot.loc[valid_seq_mask, "Tri_context"] = seq_df_not_uniprot.loc[valid_seq_mask, "Seq_dna"].apply(
        lambda x: ",".join(per_site_trinucleotide_context(x, no_flanks=True)))

    # Prepare final output
    seq_df = pd.concat((seq_df_not_uniprot, seq_df_nomane_tr)).reset_index(drop=True)
    pre_drop = len(seq_df)
    seq_df = drop_gene_duplicates(seq_df)
    logger.debug(
        "Duplicate gene removal: %s removed (from %s to %s).",
        pre_drop - len(seq_df),
        pre_drop,
        len(seq_df),
    )
    report_df = seq_df.Reference_info.value_counts().reset_index()
    report_df = report_df.rename(columns={"index" : "Source"})
    report_df.Source = report_df.Source.map({1 : "Proteins API",
                                             0 : "MANE + Ensembl API",
                                            -1 : "Backtranseq API"})
    logger_report = ", ".join([f"{v}: {c}" for (v, c) in zip(report_df.Source,
                                                              report_df.Reference_info)])
    logger.debug(f"Built of sequence dataframe completed. Retrieved {len(seq_df)} structures ({logger_report})")

    return seq_df


def get_seq_df(datasets_dir,
               output_seq_df,
               organism = "Homo sapiens",
               mane=False,
               mane_only=False,
               custom_mane_metadata_path=None,
               num_cores=1,
               mane_version=1.4):
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

    # Load Uniprot ID to HUGO and MANE to AF mapping
    pdb_dir = os.path.join(datasets_dir, "pdb_structures")
    uniprot_ids = os.listdir(pdb_dir)
    uniprot_ids = [uni_id.split("-")[1] for uni_id in list(set(uniprot_ids)) if ".pdb" in uni_id]
    logger.debug("Retrieving Uniprot ID to HUGO symbol mapping information..")
    
    if mane:
        mane_mapping, mane_mapping_not_af = get_mane_to_af_mapping(
            datasets_dir,
            uniprot_ids,
            include_not_af=True,
            mane_version=mane_version,
            custom_mane_metadata_path=custom_mane_metadata_path,
            cores=num_cores
            )
        
        uniprot_to_gene_dict = dict(zip(mane_mapping["Uniprot_ID"], mane_mapping["Gene"]))
        if custom_mane_metadata_path is not None:
            custom_symbol_map = load_custom_symbol_map(custom_mane_metadata_path)
            if custom_symbol_map:
                filled = 0
                for ens_id, symbol in custom_symbol_map.items():
                    if ens_id in uniprot_ids and (ens_id not in uniprot_to_gene_dict or pd.isna(uniprot_to_gene_dict[ens_id])):
                        uniprot_to_gene_dict[ens_id] = symbol
                        filled += 1
                if filled > 0:
                    logger.debug(
                        "Filled %s gene symbols from custom MANE samplesheet for ENSP-only entries.",
                        filled,
                    )
        missing_uni_ids = list(set(uniprot_ids) - set(mane_mapping.Uniprot_ID))
        uniprot_to_gene_dict = uniprot_to_gene_dict | uniprot_to_hugo(missing_uni_ids)
    else:
        uniprot_to_gene_dict = uniprot_to_hugo(uniprot_ids)
    
    # ---
    # # Workaround if the direct request to UniprotKB stops working (it has happened temporarily)
    # if all(pd.isna(k) for k in uniprot_to_gene_dict.keys()):
    #     logger.warning(f"Failed to retrieve Uniprot ID to HUGO symbol mapping directly from UniprotKB.")
    #     logger.warning(f"Retrying using Unipressed API client (only first HUGO symbol entry will be mapped)..")
    #     uniprot_to_gene_dict = uniprot_to_hugo_pressed(uniprot_ids)
    # ---
    
    # Create a dataframe with protein sequences
    logger.debug("Initializing sequence df..")
    seq_df = initialize_seq_df(pdb_dir, uniprot_to_gene_dict)

    if mane:
        seq_df = process_seq_df_mane(seq_df,
                                    datasets_dir,
                                    uniprot_to_gene_dict,
                                    mane_mapping, 
                                    mane_mapping_not_af,
                                    mane_only,
                                    num_cores)
    else:
        seq_df = process_seq_df(seq_df,
                                datasets_dir,
                                organism,
                                uniprot_to_gene_dict)

    # Save
    seq_df_cols = ['Gene', 'HGNC_ID', 'Ens_Gene_ID',
                   'Ens_Transcr_ID', 'Uniprot_ID', 'F',
                   'Seq', 'Chr', 'Reverse_strand', 'Exons_coord',
                   'Seq_dna', 'Tri_context','Reference_info']
    seq_df = seq_df[[col for col in seq_df_cols if col in seq_df.columns]]
    seq_df.to_csv(output_seq_df, index=False, sep="\t")
    logger.debug(f"Sequences dataframe saved in: {output_seq_df}")

    return seq_df


if __name__ == "__main__":
    raise SystemExit(
        "This module is intended to be used via the CLI: `oncodrive3d build-datasets`."
    )
