"""
Module to generate datasets necessary to run Oncodrive3D.

The build is a pipeline that perform the following tasks:
    - Download the PDB structures of the selected proteome 
      predicted by AlphaFold 2 from AlphaFold DB.
    - Merge the overlapping structures processed as fragments.
    - Extract AlphaFold model confidence (pLDDT).
    - Generate a dataframe including Uniprot_ID, Hugo Symbol,
      protein and DNA sequence. 
    - Download AlphaFold predicted aligned error (PAE) from 
      AlphaFold DB and convert the files into npy format.
    - Use the PDB structure and PAE to create maps of 
      probability of contacts (pCMAPs) for any protein of the 
      downloaded proteome with available PAE.
    - Remove unnecessary temp files (e.g., PDB structures) if 
      not specified otherwise.
"""


### TO DO ###

# Handle logging output in .sh files

#############


import daiquiri
import subprocess

import os

from scripts import __logger_name__
from scripts.datasets.af_merge import merge_af_fragments
from scripts.datasets.model_confidence import get_confidence
from scripts.datasets.seq_for_mut_prob import get_seq_df
from scripts.datasets.parse_pae import parse_pae
from scripts.datasets.prob_contact_maps import get_prob_cmaps_mp
from scripts.datasets.get_structures import get_structures

logger = daiquiri.getLogger(__logger_name__ + ".datasets.build_datasets")


def create_or_clean_dir(path_to_dir):
    """
    Create directory or clean it upon request if it already exists.
    """    
    
    if not os.path.isdir(path_to_dir):
        os.makedirs(path_to_dir)
        
    elif len(os.listdir(path_to_dir)) > 0:
        logger.warning(f"Directory {path_to_dir} already exists and is not empty.")
        overwrite = input("Do you want to overwrite it? (y/n): ")
        if overwrite.lower() in ["y", "yes"]:
            clean_files = f"rm -rf {os.path.join(path_to_dir, '*.csv')} {os.path.join(path_to_dir, '*.json')}"
            clean_pae = ["rm", "-rf", f"{path_to_dir}/pae/"]
            clean_pdb = ["rm", "-rf", f"{path_to_dir}/pdb_structures/"]
            clean_pcmaps = ["rm", "-rf", f"{path_to_dir}/prob_cmaps/"]
            subprocess.run(clean_files, shell=True)
            subprocess.run(clean_pae) 
            subprocess.run(clean_pdb) 
            subprocess.run(clean_pcmaps) 
            logger.info(f"Dataset files in {path_to_dir} have been removed.")
        else:
            logger.warning(f"Dataset files have not been removed and you will be asked to confirm overwrite for each file.")
    

def build(output_datasets, 
          organism, 
          uniprot_to_hugo, 
          num_cores, 
          af_version, 
          keep_pdb_files, 
          verbose):
    """
    Build datasets necessary to run Oncodrive3D.
    """
    
    # Paths
    dir_path = os.path.abspath(os.path.dirname(__file__))
    output_datasets = output_datasets if not None else f"{dir_path}/../../datasets"
    create_or_clean_dir(output_datasets)
        
    # Download PDB structures
    logger.info("Downloading AF predicted structures..")
    get_structures(f"{output_datasets}/pdb_structures", 
                    organism,
                    str(af_version), num_cores)

    logger.info("Download of structures completed")
        
    # Merge fragmented structures
    logger.info("Merging fragmented structures..")
    merge_af_fragments(input_dir = f"{output_datasets}/pdb_structures")
    logger.info("Merge of structures completed")
    
    # Get model confidence
    logger.info("Extracting AF model confidence..")                          # Decide what to do with default path
    get_confidence(input = f"{output_datasets}/pdb_structures", 
                   output = f"{output_datasets}/confidence.csv")
    logger.info("Extraction of model confidence completed")
    
    # Create df including genes and proteins sequences & Hugo to Uniprot_ID mapping 
    logger.info("Generating dataframe for genes and proteins sequences..")
    get_seq_df(input_dir = f"{output_datasets}/pdb_structures", 
               output_seq_df = f"{output_datasets}/seq_for_mut_prob.csv", 
               uniprot_to_gene_dict = uniprot_to_hugo, 
               organism = organism)
    logger.info("Generation of sequences dataframe completed!")
    
    # Get PAE
    logger.info("Downloading AF predicted aligned error (PAE)..")
    get_pae = [f"{dir_path}/get_pae.sh", 
               f"{output_datasets}/pdb_structures",
               f"{output_datasets}/pae",
               str(af_version),
               str(verbose)]
    subprocess.run(get_pae, check=True)
    logger.info("Download of PAE completed!")
    
    # Parse PAE
    logger.info("Parsing PAE..")                 # Might want to add multiprocessing
    parse_pae(input = f"{output_datasets}/pae")
    logger.info("Parsing PAE completed")
    
    # Get pCAMPs
    logger.info("Generating contact probability maps (pCMAPs)..")
    get_prob_cmaps_mp(input_pdb = f"{output_datasets}/pdb_structures",
                      input_pae = f"{output_datasets}/pae",
                      output = f"{output_datasets}/prob_cmaps",
                      distance = 10,
                      num_cores = num_cores)
    logger.info("Generation pCMAPs completed")
    
    # Clean datasets
    logger.info("Cleaning datasets..")
    if not keep_pdb_files:
        clean_pdb = ["rm", "-rf", f"{output_datasets}/pdb_structures/"]
        subprocess.run(clean_pdb)
    clean_pae = ["rm", "-rf", f"{output_datasets}/pae/*.json"]
    subprocess.run(clean_pae) 
    logger.info("Datasets cleaning completed")
        
    logger.info("Datasets have been successfully built and are ready for analysis!")