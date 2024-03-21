import logging
import os

import daiquiri

from scripts import __logger_name__
from scripts.plotting.utils import clean_annot_dir
from scripts.plotting.utils import get_species
from scripts.plotting.stability_change import download_stability_change, parse_ddg_rasp
from scripts.plotting.pdb_tool import run_pdb_tool, parse_pdb_tool
from scripts.plotting.pfam import get_pfam
from scripts.plotting.uniprot_feat import get_uniprot_feat

logger = daiquiri.getLogger(__logger_name__ + ".annotations.build")

logging.getLogger('urllib3.connectionpool').setLevel(logging.WARNING)



# TODO: fix multiprocessing on DDG
# TODO: multiprocessing on pdb_tool

# =================
# BUILD ANNOTATIONS 
# ================= 

def get_annotations(data_dir,
                    output_dir,
                    #path_pdb_tool_sif,
                    organism,
                    cores):
    """
    Main function to build annotations to generate annotated plots.
    """

    # Empty directory
    clean_annot_dir(output_dir, 'd')

    # Download DDG
    species = get_species(organism)
    logger.info(f"Downloading stability change...")
    if species == "Homo sapiens":
        ddg_output = os.path.join(output_dir, "stability_change")
        if not os.path.isdir(ddg_output):
            os.makedirs(ddg_output)
            logger.debug(f'mkdir {ddg_output}')
        file_path = download_stability_change(ddg_output, cores)
        logger.info(f"Download completed!")
        
        ## TODO: Optimize DDG parsing 
        ##       - only one protein is allocated to one process every time
        ##       - a list of proteins should be allocated instead
        
        # Parsing DDG
        logger.info(f"Parsing stability change...")
        parse_ddg_rasp(file_path, ddg_output, cores)
        logger.info(f"Parsing completed!")
    else:
        logger.warning(f"Currently, stability change annotation is not available for {species} but only for Homo sapiens: Skipping...")
    
    ## TODO: Enable multiprocessing for PDB_Tool
    ## TODO: Evaluate the possibility of installing PDB_Tool instead of using container
    
    # Run PDB_Tool
    logger.info(f"Extracting pACC and 2° structure...")
    path_pdb_structure = os.path.join(data_dir, "pdb_structures")
    pdb_tool_output = run_pdb_tool(input_dir=path_pdb_structure, output_dir=output_dir)
    logger.info(f"Extraction completed!")
    
    # Parse PDB_Tool
    logger.info(f"Parsing pACC and 2° structures...")
    parse_pdb_tool(input_dir=pdb_tool_output, output_dir=output_dir)
    logger.info(f"Parsing completed!")
    
    # Get Pfam annotations
    logger.info(f"Downloading and parsing Pfam...")
    get_pfam(os.path.join(output_dir, "pfam.tsv"))
    logger.info(f"Completed!")
    
    # Get Uniprot features
    logger.info(f"Downloading and parsing Features from Uniprot...")
    get_uniprot_feat(input_seq_df = os.path.join(data_dir, "seq_for_mut_prob.tsv"), 
                     output_tsv = os.path.join(output_dir, "uniprot_feat.tsv"))
    logger.info(f"Completed!")
    
    