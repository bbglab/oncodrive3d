import os
import tarfile
from scripts import __logger_name__
from scripts.datasets.utils import calculate_hash
import daiquiri
from pypdl import Downloader

logger = daiquiri.getLogger(__logger_name__ + ".build_datasets")

CHECKSUM = 'bf62d5402cb1c4580d219335a9af1ac831416edfbf2892391c8197a8356091f2'

# Usage

def extract_file(file_path, path):
     
     with tarfile.open(file_path, "r") as tar:
            tar.extractall(path)
            logger.debug(f'Extracted { int(len(tar.getnames())/2)} structure')


def download_file(url: str, destination: str, threads: int) -> None:
    """
    Downloads a file from a URL and saves it to the specified destination.

    Args:
        url (str): The URL of the file to download.
        destination (str): The local path where the file will be saved.
    """
    num_connections = 40 if threads > 40 else threads

    if os.path.exists(destination):
        logger.info(f"File {destination} already exists. Skipping download.")
    else:
        logger.debug(f'download from {url}')
        dl = Downloader()
        dl.start(url, destination, num_connections=num_connections)

    logger.debug('Asserting integrity of file:')
    try:
        assert CHECKSUM == calculate_hash(destination)
        logger.debug('File integrity check: PASS')
    except Exception as e:
        logger.critical('File integrity check: FAIL')
        logger.critical('error: ', e) 

def get_structures(path: str, species: str = 'human', af_version: str = '4', threads: int = 1) -> None:
    """
    Downloads AlphaFold predicted structures for a given organism and version.

    Args:
        path (str): Path where to save PDB structures.
        species (str, optional): Species (human (default) or mouse). Defaults to 'human'.
        af_version (str, optional): AlphaFold 2 version (4 as default). Defaults to '4'.
        verbose (str, optional): Verbose (True (default) or False). Defaults to 'False'.

    Example:
        get_structures('../../datasets/pdb_structures', species='human', verbose='True')
    """

    logger.info(f"Selected species: {species}")
    logger.info(
        f"Proteome to download: {'UP000005640_9606_HUMAN_v' + af_version}")

    if not os.path.isdir(path):
        os.makedirs(path)
        logger.debug(f'mkdir {path}')

    proteome = f"UP000005640_9606_HUMAN_v{af_version}" if species == "human" else f"UP000000589_10090_MOUSE_v{af_version}"
    af_url = f"https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/{proteome}.tar"
    file_path = os.path.join(path, f"{proteome}.tar")

    try:
        ## STEP1 --- Download file
        logger.info(f'Download to {file_path}')
        download_file(af_url, file_path, threads)

        ## STEP2 --- Extract structures
        logger.info(f'Extracting {file_path}')
        extract_file(file_path, path)

       
        # os.remove(file_path)

        logger.info('Download structure: SUCCESS')
        logger.info(f"Structures downloaded in directory {path}")

    except Exception as e:
        logger.error('Download structure: FAIL')
        logger.error(f"Error while downloading structures: {e}")
        raise e
