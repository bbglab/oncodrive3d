import logging
import os

import daiquiri
import time
import subprocess
import shutil

from scripts import __logger_name__
from scripts.datasets.utils import calculate_hash, download_single_file, extract_tar_file

logger = daiquiri.getLogger(__logger_name__ + ".build.AF-pdb")

logging.getLogger('urllib3.connectionpool').setLevel(logging.WARNING)

CHECKSUM = {
    "UP000005640_9606_HUMAN_v4": "bf62d5402cb1c4580d219335a9af1ac831416edfbf2892391c8197a8356091f2",
    "UP000000589_10090_MOUSE_v4" : "eb6c529c8757d511b75f4856c1a789378478e6255a442517ad8579708787bbab",
    "mane_overlap_v4" : "c01e9b858c5415cfe2eae7e52a561aa8a872ba0d5d4891ba0cec327b3af49d69"
}


def assert_proteome_integrity(file_path, proteome):

    if proteome in CHECKSUM.keys():
        logger.debug('Asserting integrity of file...')
        try:
            assert CHECKSUM[proteome] == calculate_hash(file_path)
            logger.debug('File integrity check: PASS')
            return "PASS"
        except Exception as e:
            logger.critical('File integrity check: FAIL')
            logger.critical(f'error: {e}') 
            return "FAIL"
    else:
        logger.warning("Assertion skipped: Proteome checksum not in records.")
        return "PASS"
    

def mv_mane_pdb(path_datasets, pdb_dir, mane_pdb_dir) -> None:
    """
    Move AF structures with overlap with MANE Select 
    transcripts to directory with AF structures from 
    human proteome. Overwrite any overlapping AF ID.
    """
        
    path_pdb = os.path.join(path_datasets, pdb_dir)
    path_mane_pdb = os.path.join(path_datasets, mane_pdb_dir)
    if not os.path.exists(path_mane_pdb):
        os.makedirs(path_mane_pdb)
    
    # Move MANE structures
    for filename in [file for file in os.listdir(path_mane_pdb) if file.endswith(".pdb.gz") or file.endswith(".pdb")]:
        pdb = os.path.join(path_pdb, filename)
        pdb_mane = os.path.join(path_mane_pdb, filename)
        shutil.move(pdb_mane, pdb)
        
    # Move MANE metadata files
    for filename in [file for file in os.listdir(os.path.join(path_mane_pdb)) if file.endswith(".csv") or file.endswith("readme.txt")]:
        source_file = os.path.join(path_mane_pdb, filename)
        dest_file = os.path.join(path_datasets, f"mane_{filename}")
        shutil.move(source_file, dest_file)
        

def get_structures(path: str, 
                   species: str = 'Homo sapiens', 
                   mane: bool = False,
                   af_version: str = '4', 
                   threads: int = 1, 
                   max_attempts: int = 15) -> None:
    """
    Downloads AlphaFold predicted structures for a given organism and version.

    Args:
        path (str): Path where to save PDB structures.
        species (str, optional): Species (human (default) or mouse). Defaults to 'human'.
        af_version (str, optional): AlphaFold 2 version (4 as default). Defaults to '4'.
        verbose (str, optional): Verbose (True (default) or False). Defaults to 'False'.

    Example:
        get_structures('datasets/pdb_structures', species='human', verbose='True')
    """

    logger.info(f"Selected species: {species}")

    if not os.path.isdir(path):
        os.makedirs(path)
        logger.debug(f'mkdir {path}')
    
    # Select proteome
    if mane:
        if species == "Homo sapiens":
            proteome = f"mane_overlap_v{af_version}"
        else:
            raise RuntimeError(f"Structures with MANE transcripts overlap are available only for 'Homo sapiens'. Exiting...")
    else:
        if species == "Homo sapiens":
            proteome = f"UP000005640_9606_HUMAN_v{af_version}"
        elif species == "Mus musculus": 
            proteome = f"UP000000589_10090_MOUSE_v{af_version}"
        else:
            raise RuntimeError(f"Failed to recognize '{species}' as organism. Currently accepted ones are 'Homo sapiens' and 'Mus musculus'. Exiting...")
            
    logger.debug(f"Proteome to download: {proteome}")
    af_url = f"https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/{proteome}.tar"
    file_path = os.path.join(path, f"{proteome}.tar")

    try:
        ## STEP1 --- Download file
        attempts = 0
        status = "INIT"
        while status != "PASS":
            download_single_file(af_url, file_path, threads)
            status = assert_proteome_integrity(file_path, proteome)
            attempts += 1
            if attempts >= max_attempts:
                raise RuntimeError(f"Failed to download with integrity after {max_attempts} attempts. Exiting...")
            time.sleep(10)
        
        ## STEP2 --- Extract structures
        logger.info(f'Extracting {file_path}')
        extract_tar_file(file_path, path)

        logger.info('Download structure: SUCCESS')
        logger.debug(f"Structures downloaded in directory {path}")

    except Exception as e:
        logger.error('Download structure: FAIL')
        logger.error(f"Error while downloading structures: {e}")
        raise e
