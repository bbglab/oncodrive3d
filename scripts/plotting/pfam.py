import logging
import os

import daiquiri
import pandas as pd
import subprocess

from scripts import __logger_name__
logger = daiquiri.getLogger(__logger_name__ + ".annotations.stability_change")

logging.getLogger('urllib3.connectionpool').setLevel(logging.WARNING)



def get_pfam(output_path):
    """
    Download and parse Pfam coordinates, name, description, 
    and Pfam ID to Transcript ID mapping.
    """
    
    try:
        # Pfam coordinates
        logger.debug("Downloading and parsing Pfam coordinates...")
        url_query = 'http://www.ensembl.org/biomart/martservice?query='
        query = '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_transcript_id" /><Attribute name = "pfam_start" /><Attribute name = "pfam_end" /><Attribute name = "pfam" /></Dataset></Query>'
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
        pfam.to_csv(output_path, index=False)
        
        # Delete temp files
        os.remove("pfam_coordinates.tsv")
        os.remove("pfam_id.tsv.gz")
        
    except Exception as e:
        logger.error('Download Pfam: FAIL')
        logger.error(f"Error while downloading Pfam: {e}")