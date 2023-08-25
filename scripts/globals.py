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
        log_dir = os.path.join(click.get_current_context().params['output_path'], 'log/')
        command_name = click.get_current_context().command.name

        if command_name == 'run':
            fname = f'{click.get_current_context().params["cohort"]}_{DATE}.log'
        else: 
            fname = f"{command_name}_{DATE}.log"

        os.makedirs(log_dir, exist_ok=True)
        
        level = logging.DEBUG if click.get_current_context().params['verbose'] else logging.INFO

        formatter = daiquiri.formatter.ColorFormatter(fmt=FORMAT)
        
        daiquiri.setup(level=level, outputs=(
            daiquiri.output.Stream(formatter=formatter), 
            daiquiri.output.File(filename=os.path.join(log_dir, fname), formatter=formatter)
        ))
        
        logger.debug(f'Log path: {os.path.join(str(log_dir), fname)}')
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