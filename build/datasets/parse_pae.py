import numpy as np
import os
import json
import click
import re
from progressbar import progressbar


###################################### EXAMPLE USAGE ################################################

# python3 parse_pae.py -i /workspace/projects/clustering_3d/clustering_3d/datasets_frag/pae

#####################################################################################################

##### NOTE #####
# 1) In the final pipeline the PAE files in original format must be deleted after parsing
# 2) I might want to enable multiprocessing in this step
################


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


@click.command(context_settings=dict(help_option_names=['-h', '--help']),
               help='Convert all predicted aligned error from .json dict to .npy array.')
@click.option("-i", "--input_dir", type=click.Path(exists=True), required=True, 
              help="Path to directory including the json files (predicted aligned error)")
@click.option("-o", "--output", help="Path to output directory", type=str)
def main(input_dir, output):
    
    if output is None:
        output = input_dir
        
    # Create necessary folder
    if not os.path.exists(output):
        os.makedirs(output)
        
    path_files = get_pae_path_list_from_dir(input_dir)

    for path in progressbar(path_files):
        np.save(path.replace(".json", ".npy"), json_to_npy(path))


if __name__ == "__main__":
    main()