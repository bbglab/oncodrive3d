import logging
import os

import daiquiri

from scripts import __logger_name__
from scripts.plotting.utils import clean_annot_dir
from scripts.plotting.stability_change import download_stability_change, parse_ddg_rasp
from scripts.plotting.pdb_tool import run_pdb_tool, parse_pdb_tool
from scripts.plotting.pfam import get_pfam

logger = daiquiri.getLogger(__logger_name__ + ".annotations.build")

logging.getLogger('urllib3.connectionpool').setLevel(logging.WARNING)



# TODO: fix multiprocessing on DDG
# TODO: multiprocessing on pdb_tool

# =================
# BUILD ANNOTATIONS 
# ================= 

def get_annotations(data_dir,
                    path_pdb_tool_sif,
                    output_dir,
                    cores,
                    verbose):
    """
    Main function to build annotations to generate annotated plots.
    """

    # Empty directory
    clean_annot_dir(output_dir, 'd')

    # Download DDG
    logger.info(f"Downloading stability change...")
    ddg_output = f"{output_dir}/stability_change"
    if not os.path.isdir(ddg_output):
        os.makedirs(ddg_output)
        logger.debug(f'mkdir {ddg_output}')
    file_path = download_stability_change(ddg_output, cores)
    logger.info(f"Download completed!")
    
    # Parsing DDG
    logger.info(f"Parsing stability change...")
    parse_ddg_rasp(file_path, ddg_output, cores)
    logger.info(f"Parsing completed!")
    
    ## TODO: Enable multiprocessing for PDB_Tool
    ## TODO: Find a way to avoid using container for PDB_Tool?
    
    # Run PDB_Tool
    logger.info(f"Extracting pACC and 2° structure...")
    path_pdb_structure = f"{data_dir}/pdb_structures"
    pdb_tool_output = run_pdb_tool(path_pdb_tool_sif, inpur_dir=path_pdb_structure, output_dir=output_dir)
    logger.info(f"Extraction completed!")
    
    # Parse PDB_Tool
    logger.info(f"Parsing pACC and 2° structures...")
    parse_pdb_tool(input_dir=pdb_tool_output, output_dir=output_dir)
    logger.info(f"Parsing completed!")
    
    # Get Pfam annotations
    logger.info(f"Downloading and parsing Pfam...")
    get_pfam(f"{output_dir}/pfam.tsv")
    logger.info(f"Completed!")