"""
Simple module to convert the predicted aligned error (PAE) 
files produced by AlphaFold2 from json to npy format.
"""


import numpy as np
import os
import json
import re
from progressbar import progressbar


def get_pae_path_list_from_dir(path_dir):
    """
    Takes as input the path of a directory and it 
    outputs a list of paths of the contained PAE files.
    """

    pae_files = os.listdir(path_dir)
    pae_files = [f"{path_dir}/{f}" for f in pae_files if re.search('-predicted_aligned_error.json$', f) is not None]
    
    return pae_files


def json_to_npy(path):
    
    with open(path) as f:
        pae = json.load(f)
        
    return np.array(pae[0]['predicted_aligned_error'])


def parse_pae(input, output=None):
    
    if output is None:
        output = input
        
    # Create necessary folder
    if not os.path.exists(output):
        os.makedirs(output)
        
    path_files = get_pae_path_list_from_dir(input)

    for path in progressbar(path_files):
        np.save(path.replace(".json", ".npy"), json_to_npy(path))