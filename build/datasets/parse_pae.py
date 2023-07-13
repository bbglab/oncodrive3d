import numpy as np
import os
import json
import argparse
import re
from progressbar import progressbar


###################################### EXAMPLE USAGE ################################################

# python3 parse_pae.py -i /workspace/projects/clustering_3d/clustering_3d/datasets_frag/pae

#####################################################################################################

##### NOTE #####
# I might want to add a final check that ensure that there are no truncated structures
################

def init_parser():
    """
    Initialize parser for the main function.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Path to directory including the json files (predicted aligned error)", type=str, required=True) 
    parser.add_argument("-o", "--output", help="Path to the output dir", type=str) 
    return parser.parse_args()


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


def main():
    
    args = init_parser()
    input = args.input
    output = args.output
    
    if output is None:
        output = input
        
    # Create necessary folder
    if not os.path.exists(output):
        os.makedirs(output)
        
    path_files = get_pae_path_list_from_dir(input)

    for path in progressbar(path_files):
        np.save(path.replace(".json", ".npy"), json_to_npy(path))

if __name__ == "__main__":
    main()