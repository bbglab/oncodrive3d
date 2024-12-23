import logging
import os
import subprocess

import daiquiri
import pandas as pd

from scripts import __logger_name__

logger = daiquiri.getLogger(__logger_name__ + ".plotting.pfam")

logging.getLogger('urllib3.connectionpool').setLevel(logging.WARNING)



def add_pfam_metadata(pfam, seq_df):
    """
    Add Ensembl transcript and gene info and rename cols to
    be merged with Uniprot features dataframe.
    """

    # Add metadata to PFAM
    pfam = seq_df[["Gene", "Uniprot_ID", "Ens_Transcr_ID", "Ens_Gene_ID"]].merge(
        pfam, how="left", on=["Ens_Transcr_ID", "Ens_Gene_ID"])
    pfam = pfam.dropna(how="all", subset=["Pfam_start", "Pfam_end"]).reset_index(drop=True)

    # Prepare to merge
    pfam["Type"] = "DOMAIN"
    pfam["Evidence"] = "Pfam"
    pfam = pfam.rename(columns={"Pfam_start" : "Begin",
                                "Pfam_end" : "End",
                                "Pfam_name" : "Description",
                                "Pfam_description" : "Full_description"})

    return pfam


def get_pfam(seq_df, output_tsv, organism):
    """
    Download and parse Pfam coordinates, name, description,
    and Pfam ID to Transcript ID mapping.
    """

    status = "INIT"
    i = 0
    if organism == "Homo sapiens":
        ensembl_gene_dataset = "hsapiens_gene_ensembl"
    elif organism == "Mus musculus":
        ensembl_gene_dataset = "mmusculus_gene_ensembl"
    else:
        logger.error(f"Invalid organism: {organism}. Expected 'Homo sapiens' or 'Mus musculus'.")
        raise ValueError(f"Invalid organism: {organism}. Must be 'Homo sapiens' or 'Mus musculus'.")

    while status != "PASS":
        if i < 5:
            try:
                # Pfam coordinates
                logger.debug("Downloading and parsing Pfam coordinates...")
                url_query = 'http://jan2024.archive.ensembl.org/biomart/martservice?query='
                query = f'<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "{ensembl_gene_dataset}" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_transcript_id" /><Attribute name = "pfam_start" /><Attribute name = "pfam_end" /><Attribute name = "pfam" /></Dataset></Query>'
                url = url_query + query
                command = [f"wget", "-q", "-O", "pfam_coordinates.tsv", url]
                subprocess.run(command)
                pfam = pd.read_csv("pfam_coordinates.tsv", sep="\t", header=None)
                pfam.columns = ["Ens_Gene_ID", "Ens_Transcr_ID", "Pfam_start", "Pfam_end", "Pfam_ID"]

                # ID database
                logger.debug("Downloading and parsing Pfam ID database...")
                url = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/database_files/pfamA.txt.gz"
                command = ["wget", "-q", "-O", "pfam_id.tsv.gz", url]
                subprocess.run(command)
                pfam_id = pd.read_csv("pfam_id.tsv.gz", compression='gzip', sep='\t', header=None).iloc[:,[0,1,3]]
                pfam_id.columns = "Pfam_ID", "Pfam_name", "Pfam_description"

                # Merge and save
                pfam = pfam.merge(pfam_id, how="left", on="Pfam_ID")
                pfam = pfam.dropna(how="all", subset=["Pfam_start", "Pfam_end"]).reset_index(drop=True)
                pfam = add_pfam_metadata(pfam, seq_df)
                pfam.to_csv(output_tsv, index=False, sep="\t")

                # Delete temp files
                os.remove("pfam_coordinates.tsv")
                os.remove("pfam_id.tsv.gz")
                status = "PASS"

                return pfam

            except Exception as e:
                status = "FAIL"
                logger.warning(f"Error while downloading Pfam: {e}")
                logger.warning("Retrying download...")
                i += 1
        else:
            logger.error(f'Download Pfam: {status}')