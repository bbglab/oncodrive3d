import logging
import os

import daiquiri
import subprocess
import click

from scripts import __logger_name__

logger = daiquiri.getLogger(__logger_name__ + ".annotations.utils")

logging.getLogger('urllib3.connectionpool').setLevel(logging.WARNING)



# =====
# Utils
# =====


def clean_annotations_dir(path: str, loc: str) -> None:
    """
    Clean the annotations directory by removing specific files 
    and subdirectories.

    Args:
        path (str): Path to the directory to be cleaned.
    """

    if loc == "d":

        clean_files = f"rm -rf {os.path.join(path, '*.csv')} {os.path.join(path, '*.json')} {os.path.join(path, '.*.txt')}"
        clean_ddg = ["rm", "-rf", os.path.join(path, "stability_change")]
        clean_pdbtool = ["rm", "-rf", os.path.join(path, "pdb_tool")]
        #clean_log = ["rm", "-rf", os.path.join(path, "log")]

        logger.debug(clean_files)
        subprocess.run(clean_files, shell=True)

        logger.debug(' '.join(clean_ddg))
        subprocess.run(clean_ddg)
        
        logger.debug(' '.join(clean_pdbtool))
        subprocess.run(clean_pdbtool)

        # logger.debug(' '.join(clean_log))
        # subprocess.run(clean_log)

    elif loc == "r":
        # TODO: implement cleaning function for output
        pass


def clean_annot_dir(path: str, loc: str = 'd') -> None:
    """
    Clean it upon request if it already exists.

    Args:
        path (str): Path to the directory to be created or cleaned.
    """

    if os.listdir(path) != ['log']:
        logger.warning(f"Directory {path} already exists and is not empty.")

        overwrite = "y" if click.get_current_context().params['yes'] else input("Clean existing directory? (y/n): ")
        while overwrite.lower() not in ["y", "yes", "n", "no"]:
            print("Please choose yes or no")
            overwrite = input("Clean existing directory? (y/n): ")

        if overwrite.lower() in ["y", "yes"]:
            clean_annotations_dir(path, loc)
            logger.info(f"Dataset files in {path} have been removed.")
        else:
            logger.warning(f"Dataset files in {path} have not been removed.")
    else:
        pass