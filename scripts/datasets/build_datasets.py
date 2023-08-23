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

# Convert print into logs in the other modules
# How the logger from different modules understand the level?
# > Check that the level is defined correctly across modules
# > Remove verbose inside functions if it is already defined in the level
# Handle warnings
# Handle cases where the datasets folder already exist
# > Overwrite? Maybe ask? If so, remove everything already present
# Handle logging output in .sh files

#############

from scripts import __logger_name__

import logging
import subprocess
import os

from scripts.datasets.af_merge import merge_af_fragments
from scripts.datasets.model_confidence import get_confidence
from scripts.datasets.seq_for_mut_prob import get_seq_df
from scripts.datasets.parse_pae import parse_pae
from scripts.datasets.prob_contact_maps import get_prob_cmaps_mp

logger = logging.getLogger(__logger_name__ + ".build_datasets")



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
    if not os.path.isdir(output_datasets):
        os.makedirs(output_datasets)

    # Download PDB structures
    logger.info("Downloading AF predicted structures..")
    download_pdb = [f"{dir_path}/get_structures.sh", 
                    f"{output_datasets}/pdb_structures", 
                    organism,
                    str(af_version), 
                    str(verbose)]
    subprocess.run(download_pdb, check=True)                                                     ### >> Is this the correct way of handling error? <<
    logger.info("Download of structures completed", flush=True)
        
    # Merge fragmented structures
    logger.info("Merging fragmented structures..", flush=True)
    merge_af_fragments(input_dir = f"{output_datasets}/pdb_structures")
    logger.info("Merge of structures completed", flush=True)
    
    # Get model confidence
    logger.info("Extracting AF model confidence..", flush=True)                                    # Decide what to do with default path
    get_confidence(input = f"{output_datasets}/pdb_structures", 
                   output = f"{output_datasets}/confidence.csv",
                   verbose = verbose)
    logger.info("Extraction of model confidence completed", flush=True)
    
    # Create df including genes and proteins sequences & Hugo to Uniprot_ID mapping 
    logger.info("Generating dataframe for genes and proteins sequences..", flush=True)
    get_seq_df(input_dir = f"{output_datasets}/pdb_structures", 
               output_seq_df = f"{output_datasets}/seq_for_mut_prob.csv", 
               uniprot_to_gene_dict = uniprot_to_hugo, 
               organism = organism,
               verbose = verbose)
    logger.info("Generation of sequences dataframe completed", flush=True)
    
    # Get PAE
    logger.info("Downloading AF predicted aligned error (PAE)..", flush=True)
    get_pae = [f"{dir_path}/get_pae.sh", 
               f"{output_datasets}/pdb_structures",
               f"{output_datasets}/pae",
               str(af_version),
               str(verbose)]
    subprocess.run(get_pae, check=True)
    logger.info("Download of PAE completed", flush=True)
    
    # Parse PAE
    logger.info("Parsing PAE..", flush=True)                 # Might want to delete original PAE & might want to add multiprocessing
    parse_pae(input = f"{output_datasets}/pae")
    logger.info("Parsing PAE completed", flush=True)
    
    # Get pCAMPs
    logger.info("Generating contact probability maps (pCMAPs)..", flush=True)
    get_prob_cmaps_mp(input_pdb = f"{output_datasets}/pdb_structures",
                      input_pae = f"{output_datasets}/pae",
                      output = f"{output_datasets}/prob_cmaps",
                      distance = 10,
                      num_cores = num_cores,
                      verbose = verbose)
    logger.info("Generation pCMAPs completed", flush=True)
    
    # Clean datasets
    logger.info("Cleaning datasets..", flush=True)
    if not keep_pdb_files:
        clean_pdb = ["rm", "-rf", f"{output_datasets}/pdb_structures/"]
        subprocess.run(clean_pdb)
    clean_pae = ["rm", "-rf", f"{output_datasets}/pae/*.json"]
    subprocess.run(clean_pae) 
    logger.info("Datasets cleaning completed", flush=True)
        
    logger.info("Building datasets completed", flush=True)