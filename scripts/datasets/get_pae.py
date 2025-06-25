import os

import daiquiri
import requests
import time
from typing import Optional
import concurrent.futures
from tqdm import tqdm

from scripts import __logger_name__

logger = daiquiri.getLogger(__logger_name__ + ".build.PAE")


def download_pae(
    uniprot_id: str, 
    af_version: int, 
    output_dir: str,
    max_retries: int = 100
    ) -> None:
    """
    Download Predicted Aligned Error (PAE) file from AlphaFold DB.

    Args:
        uniprot_id: Uniprot ID of the structure.
        af_version: AlphaFold 2 version.
        output_dir: Output directory where to download the PAE files.
        max_retries: Break the loop if the download fails too many times.
    """

    file_path = os.path.join(output_dir, f"{uniprot_id}-F1-predicted_aligned_error.json")
    download_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-predicted_aligned_error_v{af_version}.json"

    i = 0
    status = "INIT"
    while status != "FINISHED":
        i += 1
        if i > max_retries:
            logger.warning(f"Reached {max_retries} attempts; proceeding without download.")
            break
        if status != "INIT":
            time.sleep(30)
        try:
            response = requests.get(download_url, timeout=30)
            content = response.content
            if content.endswith(b'}]') and not content.endswith(b'</Error>'):
                with open(file_path, 'wb') as output_file:
                    output_file.write(content)
                status = "FINISHED"
        except requests.exceptions.RequestException as e:
            status = "ERROR"
            if i % 10 == 0:
                logger.debug(f"Request failed {e}: Retrying")


def get_pae(
    input_dir: str, 
    output_dir: str, 
    num_cores: int, 
    af_version: int = 4,
    custom_pdb_dir: Optional[str] = None
    ) -> None:
    """
    Download Predicted Aligned Error (PAE) files for all non-fragmented PDB
    structures in the input directory.

    Args:
        input_dir (str): Input directory including the PDB structures.
        output_dir (str): Output directory where to download the PAE files.
        num_cores (int): Number of cores for multithreading download.
        af_version (int): AlphaFold 2 version (default is 4).
        custom_pdb_dir (str | None): Directory including provided custom PDB structures.
    """

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    checkpoint = os.path.join(output_dir, '.checkpoint.txt')
    if os.path.exists(checkpoint):
        logger.debug("PAE already downloaded: Skipping...")
        return
    
    # Do not download PAE of custom provided structures
    if custom_pdb_dir is not None:
        custom_uniprot_ids = [fname.split('.')[0] for fname in os.listdir(custom_pdb_dir) if fname.endswith('.pdb')]

    pdb_files = [file for file in os.listdir(input_dir) if file.startswith("AF-") and file.endswith(f"-model_v{af_version}.pdb.gz")]
    uniprot_ids = [pdb_file.split("-")[1] for pdb_file in pdb_files]
    uniprot_ids = [uni_id for uni_id in uniprot_ids if uni_id not in custom_uniprot_ids]

    with concurrent.futures.ThreadPoolExecutor(max_workers=num_cores) as executor:
        tasks = [executor.submit(download_pae, uniprot_id, af_version, output_dir) for uniprot_id in uniprot_ids]

        for _ in tqdm(concurrent.futures.as_completed(tasks), total=len(tasks), desc="Downloading PAE"):
            pass

    with open(checkpoint, "w") as f:
        f.write('')

    logger.info('Download of PAE completed!')