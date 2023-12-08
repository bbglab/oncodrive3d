import logging
import os

import daiquiri

from scripts import __logger_name__
from scripts.plotting.utils import clean_annot_dir
from scripts.plotting.stability_change import download_stability_change, parse_ddg_rasp
from scripts.plotting.pdb_tool import run_pdb_tool, load_all_pdb_tool_files, pdb_tool_to_3s_sse

logger = daiquiri.getLogger(__logger_name__ + ".annotations.build")

logging.getLogger('urllib3.connectionpool').setLevel(logging.WARNING)



# TODO: fix multiprocessing on DDG
# TODO: multiprocessing on pdb_tool

# =================
# BUILD ANNOTATIONS 
# ================= 

def get_annotations(path_pdb_structure,
                    path_pdb_tool_sif,
                    output_dir,
                    cores,
                    verbose):
    """
    Main function to build annotations to generate annotated plots.
    """
    
    logger.info(f"Path to PDB structures: {path_pdb_structure}")
    logger.info(f"Output to PDB_Tool SIF: {path_pdb_tool_sif}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Cores: {cores}")
    logger.info(f"Verbose: {bool(verbose)}")
    logger.info(f'Log path: {os.path.join(output_dir, "log")}')
    logger.info("")

    # Empty directory
    clean_annot_dir(output_dir, 'd')

    # Download DDG
    logger.info(f"Downloading stability change...")
    ddg_output = f"{output_dir}/stability_change"
    if not os.path.isdir(ddg_output):
        os.makedirs(ddg_output)
        logger.debug(f'mkdir {ddg_output}')
    file_path = download_stability_change(ddg_output, cores)
    logger.info(f"Download stability change completed!")
    
    # Parsing DDG
    logger.info(f"Parsing stability change...")
    parse_ddg_rasp(file_path, ddg_output, cores)
    logger.info(f"Parsing stability change completed!")
    
    # # TODO: Enable multiprocessing for PDB_Tool 
    
    # Run PDB_Tool
    logger.info(f"Extracting pACC and 2° structure...")
    pdb_tool_output = f"{output_dir}/pdb_tool"
    if not os.path.isdir(pdb_tool_output):
        os.makedirs(pdb_tool_output)
        logger.debug(f'mkdir {pdb_tool_output}')
    run_pdb_tool(path_pdb_tool_sif, i=path_pdb_structure, o=pdb_tool_output)
    logger.info(f"Extraction completed!")
    
    # Parse PDB_Tool
    logger.info(f"Parsing pACC and 2° structures...")
    pdb_tool_df = load_all_pdb_tool_files(pdb_tool_output)
    logger.debug(f"Deleting {pdb_tool_output}")
    os.rmdir(pdb_tool_output)
    pdb_tool_df = pdb_tool_to_3s_sse(pdb_tool_df)
    pdb_tool_df = pdb_tool_df.drop(columns=["CLE", "ACC", "CNa", "CNb"])
    pdb_tool_df.to_csv(f"{output_dir}/pdb_tool_df.csv", index=False)
    logger.info(f"Parsing completed!")