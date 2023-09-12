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
# Suppress warnings in merge_af.py for atoms

#############


import os

import daiquiri

from scripts import __logger_name__
from scripts.datasets.af_merge import merge_af_fragments
from scripts.datasets.get_pae import get_pae
from scripts.datasets.get_structures import get_structures
from scripts.datasets.model_confidence import get_confidence
from scripts.datasets.parse_pae import parse_pae
from scripts.datasets.prob_contact_maps import get_prob_cmaps_mp
from scripts.datasets.seq_for_mut_prob import get_seq_df
from scripts.globals import clean_dir, clean_temp_files

logger = daiquiri.getLogger(__logger_name__ + ".build")

def build(output_datasets,
          organism,
          distance_threshold,
          uniprot_to_hugo,
          num_cores,
          af_version,
          keep_pdb_files,
          ):
    """
    Build datasets necessary to run Oncodrive3D.
    """

    # empty directory
    clean_dir(output_datasets, 'd')

    # Download PDB structures
    logger.info("Downloading AlphaFold (AF) predicted structures...")
    get_structures(path=os.path.join(output_datasets,"pdb_structures"),
                   species=organism,
                   af_version=str(af_version), 
                   threads=num_cores)

    logger.info("Download of structures completed!")

    # Merge fragmented structures
    logger.info("Merging fragmented structures...")
    merge_af_fragments(input_dir=os.path.join(output_datasets,"pdb_structures"), 
                       gzip=True)
    logger.info("Merge of structures completed!")

    # Get model confidence
    # Decide what to do with default path
    logger.info("Extracting AF model confidence...")
    get_confidence(input=os.path.join(output_datasets, "pdb_structures"),
                   output_dir=os.path.join(output_datasets))

    # Create df including genes and proteins sequences & Hugo to Uniprot_ID mapping
    logger.info("Generating dataframe for genes and proteins sequences...")
    get_seq_df(input_dir=os.path.join(output_datasets,"pdb_structures"),
               output_seq_df=os.path.join(output_datasets, "seq_for_mut_prob.csv"),
               uniprot_to_gene_dict=uniprot_to_hugo,
               organism=organism)
    logger.info("Generation of sequences dataframe completed!")

    # Get PAE
    logger.info("Downloading AF predicted aligned error (PAE)...")
    get_pae(input_dir=os.path.join(output_datasets,"pdb_structures"),
            output_dir=os.path.join(output_datasets,"pae"),
            af_version=str(af_version),
            )

    # Parse PAE
    # Might want to add multiprocessing
    logger.info("Parsing PAE...")
    parse_pae(input=os.path.join(output_datasets, 'pae'))
    logger.info("Parsing PAE completed!")

    # Get pCAMPs
    logger.info("Generating contact probability maps (pCMAPs)..")
    get_prob_cmaps_mp(input_pdb=os.path.join(output_datasets, "pdb_structures"),
                      input_pae=os.path.join(output_datasets, "pae"),
                      output=os.path.join(output_datasets,"prob_cmaps"),
                      distance=distance_threshold,
                      num_cores=num_cores)
    logger.info("Generation pCMAPs completed!")

    # Clean datasets
    logger.info("Cleaning datasets...")
    clean_temp_files(path=output_datasets,
                     keep_pdb_files=keep_pdb_files)
    logger.info("Datasets cleaning completed!")

    logger.info("Datasets have been successfully built and are ready for analysis!")
