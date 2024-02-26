import logging
import os

import daiquiri
import time
import subprocess

from scripts import __logger_name__
from scripts.datasets.utils import calculate_hash, download_single_file, extract_tar_file

logger = daiquiri.getLogger(__logger_name__ + ".build.AF-pdb")

logging.getLogger('urllib3.connectionpool').setLevel(logging.WARNING)

CHECKSUM = {
    "UP000005640_9606_HUMAN_v4": 'bf62d5402cb1c4580d219335a9af1ac831416edfbf2892391c8197a8356091f2',
}


def assert_integrity_human(file_path, proteome):

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


def get_structures(path: str, 
                   species: str = 'Homo sapiens', 
                   af_version: str = '4', 
                   threads: int = 1, 
                   max_attempts: int = 10) -> None:
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
            status = assert_integrity_human(file_path, proteome)
            attempts += 1
            if attempts >= max_attempts:
                raise RuntimeError(f"Failed to download with integrity after {max_attempts} attempts. Exiting...")
            time.sleep(10)
        
        ## STEP2 --- Extract structures
        logger.info(f'Extracting {file_path}')
        extract_tar_file(file_path, path)

        # os.remove(file_path)

        logger.info('Download structure: SUCCESS')
        logger.debug(f"Structures downloaded in directory {path}")

    except Exception as e:
        logger.error('Download structure: FAIL')
        logger.error(f"Error while downloading structures: {e}")
        raise e
