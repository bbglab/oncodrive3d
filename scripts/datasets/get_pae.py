import os

import daiquiri
import requests
from progressbar import progressbar

from scripts import __logger_name__

logger = daiquiri.getLogger(__logger_name__ + ".build.PAE")


def download_pae(uniprot_id: str, af_version: int, output_dir: str) -> None:
    """
    Download Predicted Aligned Error (PAE) files.

    Args:
        uniprot_id (str): Uniprot ID of the structure.
        af_version (int): AlphaFold 2 version.
        output_dir (str): Output directory where to download the PAE files.
    """
    file_path = os.path.join(output_dir, f"{uniprot_id}-F1-predicted_aligned_error.json")
    download_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-predicted_aligned_error_v{af_version}.json"
    
    response = requests.get(download_url)
    content = response.content

    if content.endswith(b'}]') and not content.endswith(b'</Error>'):
        with open(file_path, 'wb') as output_file:
            output_file.write(content)


def get_pae(input_dir: str, output_dir: str, af_version: int = 4) -> None:
    """
    Download Predicted Aligned Error (PAE) files for non-fragmented PDB structures.

    Args:
        input_dir (str): Input directory including the PDB structures.
        output_dir (str): Output directory where to download the PAE files.
        af_version (int): AlphaFold 2 version (default is 4).
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    checkpoint = os.path.join(output_dir, '.checkpoint.txt')
    if os.path.exists(checkpoint):
        logger.debug("PAE already downloaded. Skipping")

    else:
        for pdb_file in progressbar(os.listdir(input_dir)):
            if pdb_file.startswith("AF-") and pdb_file.endswith(f"-model_v{af_version}.pdb.gz"):
                uniprot_id = pdb_file.split("-")[1]
                download_pae(uniprot_id, af_version, output_dir)

        with open(checkpoint, "w") as f:
                    f.write('')

        logger.info('Download completed.')