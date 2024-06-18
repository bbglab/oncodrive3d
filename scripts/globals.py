import os
import logging
import daiquiri
import click
import shutil
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

        if command_name in ['run', 'plot', 'chimerax_plot', 'chimerax-plot']:
            cohort = click.get_current_context().params["cohort"]
            fname = f'{command_name}_{cohort if cohort != "None" else ""}_{DATE}.log'
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
    support_email = "stefano.pellegrini@irbbarcelona.org"
    banner_width = 70

    logger.info("#" * banner_width)
    logger.info(f"{'#' + ' ' * (banner_width - 2) + '#'}")
    logger.info(f"{'#' + f'Welcome to Oncodrive3D!'.center(banner_width - 2) + '#'}")
    logger.info(f"{'#' + ' ' * (banner_width - 2) + '#'}")
    logger.info(f"{'#' + initializing_text.center(banner_width - 2) + '#'}")
    logger.info(f"{'#' + ' ' * (banner_width - 2) + '#'}")
    logger.info(f"{'#' + f'Version: {version}'.center(banner_width - 2) + '#'}")
    logger.info(f"{'#' + f'Author: {author}'.center(banner_width - 2) + '#'}")
    logger.info(f"{'#' + f'Support: {support_email}'.center(banner_width - 2) + '#'}")
    logger.info(f"{'#' + ' ' * (banner_width - 2) + '#'}")
    logger.info("#" * banner_width)
    logger.info("")


# ===================
#  Clean and organize
# ===================


def copy_dir(source_dir, destination_dir):
    """
    Ccopy the entire directory.
    """

    logger.debug("Copying directory..")
    logger.debug(f"From {source_dir}")
    logger.debug(f"To {destination_dir}")
    
    if os.path.exists(source_dir):
        if os.path.isdir(source_dir):
            
            try:
                shutil.copytree(source_dir, destination_dir)
                logger.debug("Directory copied successfully!")
            except Exception as e:
                logger.warning(f"An error occurred ({e}): Skipping")
                
        else:
            logger.warning("Error while copying directory (source path is not a directory): Skipping")
    else:
        logger.warning("Error while copying directory (source path does not exist): Skipping")
        

def rm_dir(dir_path):             # TO DO: Probably, if the directory is not empty, it should ask for confirmation
    """
    Remove directory.
    """
    
    if os.path.isdir(dir_path):
        shutil.rmtree(dir_path, ignore_errors=True)
    
    
def rm_files(dir_path, ext=[".cif.gz", ".cif"]) -> None:
    """
    Remove any file with a given extension in a given directory.
    """
    
    for file in os.listdir(dir_path):
        file_path = os.path.join(dir_path, file)
        if any([file.endswith(ext_i) for ext_i in ext]):
            os.remove(file_path)


def clean_directory(path: str, loc: str) -> None:
    """
    Clean a directory by removing specific files and subdirectories.
    """

    if loc == "d":
        logger.debug(f"Cleaning {path}")
        rm_files(path, ext=[".tsv", ".csv", ".json", ".txt", ".txt.gz"])
        dirs = "pae", "pdb_structures", "pdb_structures_mane", "prob_cmaps"
        path_dirs = [os.path.join(path, d) for d in dirs]
        for path in path_dirs:
            rm_dir(path)

    elif loc == "r":
        # TODO: implement cleaning function for output
        pass


def clean_dir(path: str, loc: str = 'd', txt_file=False) -> None:
    """
    Clean it upon request if it already exists.
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
            
            
def clean_temp_files(path: str) -> None:
    """
    Clean temp files from dir after completing building the datasets. 
    """
    
    pdb_dir = os.path.join(path, "pdb_structures")
    rm_files(pdb_dir, ext=[".cif.gz", ".cif", ".tar"])
    rm_files(os.path.join(path, "pae"), ext=[".json"])
    rm_dir(os.path.join(path, "pdb_structures_mane"))
    rm_dir(os.path.join(pdb_dir, "fragmented_pdbs"))