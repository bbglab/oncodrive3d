import os
import logging
import daiquiri
import click
from datetime import datetime

from functools import wraps

from scripts import __logger_name__


logger = daiquiri.getLogger(__logger_name__)

DATE = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
FORMAT = "%(asctime)s - %(color)s%(levelname)-7s%(color_stop)s | %(name)s - %(color)s%(message)s%(color_stop)s"


# =========
#  Logging
# =========

def setup_logging_decorator(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        log_dir = os.path.join(click.get_current_context().params['output_dir'], 'log')
        command_name = click.get_current_context().command.name

        if command_name == 'run':
            cohort = click.get_current_context().params["cohort"]
            fname = f'{cohort if not "None" else command_name}_{DATE}.log'
        else: 
            fname = f"{command_name}_{DATE}.log"

        os.makedirs(log_dir, exist_ok=True)
        
        level = logging.DEBUG if click.get_current_context().params['verbose'] else logging.INFO

        formatter = daiquiri.formatter.ColorFormatter(fmt=FORMAT)
        
        daiquiri.setup(level=level, outputs=(
            daiquiri.output.Stream(formatter=formatter), 
            daiquiri.output.File(filename=os.path.join(log_dir, fname), formatter=formatter)
        ))
        
        return func(*args, **kwargs)

    return wrapper


def startup_message(version, initializing_text):
    
    author = "Biomedical Genomics Lab - IRB Barcelona"
    support_email = "stefano.pellegrini@irbbarcelona.com"
    banner_width = 70

    logger.info("#" * banner_width)
    logger.info(f"{'#' + ' ' * (banner_width - 2) + '#'}")
    logger.info(f"{'#' + f'Welcome to Oncodrive3D!'.center(banner_width - 2) + '#'}")
    logger.info(f"{'#' + ' ' * (banner_width - 2) + '#'}")
    logger.info(f"{'#' + initializing_text.center(banner_width - 2) + '#'}")
    logger.info(f"{'#' + f'Version: {version}'.center(banner_width - 2) + '#'}")
    logger.info(f"{'#' + f'Author: {author}'.center(banner_width - 2) + '#'}")
    logger.info(f"{'#' + f'Support: {support_email}'.center(banner_width - 2) + '#'}")
    logger.info(f"{'#' + ' ' * (banner_width - 2) + '#'}")
    logger.info("#" * banner_width)
    logger.info("")


# =========
#  Clean
# =========

import os
import subprocess
import logging

def clean_directory(path: str, loc: str) -> None:
    """
    Clean a directory by removing specific files and subdirectories.

    Args:
        path (str): Path to the directory to be cleaned.
    """

    if loc == "d":

        clean_files = f"rm -rf {os.path.join(path, '*.csv')} {os.path.join(path, '*.json')} {os.path.join(path, '.*.txt')}"
        clean_pae = ["rm", "-rf", os.path.join(path, "pae")]
        clean_pdb = ["rm", "-rf", os.path.join(path, "pdb_structures")]
        clean_pcmaps = ["rm", "-rf", os.path.join(path, "prob_cmaps")]

        logger.debug(clean_files)
        subprocess.run(clean_files, shell=True)

        logger.debug(' '.join(clean_pae))
        subprocess.run(clean_pae)

        logger.debug(' '.join(clean_pdb))
        subprocess.run(clean_pdb)

        logger.debug(' '.join(clean_pcmaps))
        subprocess.run(clean_pcmaps)

    elif loc == "r":
        # TODO: implement cleaning function for output
        pass


def clean_dir(path: str, loc: str = 'd') -> None:
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
            clean_directory(path, loc)
            logger.info(f"Dataset files in {path} have been removed.")
        else:
            logger.warning(f"Dataset files in {path} have not been removed.")
    else:
        pass
    
    
def clean_temp_files(path: str, keep_pdb_files: bool) -> None:
    """
    Clean temp files from dir after completing building the datasets. 

    Args:
        path (str): Path to build directory to be cleaned.
    """
    
    if not keep_pdb_files:
        clean_pdb = ["rm", "-rf", os.path.join(path, "pdb_structures")]
        logger.debug(' '.join(clean_pdb))
        subprocess.run(clean_pdb)
    clean_pae = ["rm", "-rf", os.path.join(path, "pae", "*.json")]
    logger.debug(' '.join(clean_pae))
    subprocess.run(clean_pae)
