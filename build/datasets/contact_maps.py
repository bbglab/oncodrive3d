"""
Compute contact maps related features from AlphaFold predicted structures.
The output can be either individual contact maps named with UniprotID_F 
(.npy or .csv), where F is the number of AF fragment. It can generate, or
a single dictionary having UniprotID_F as keys and contact maps as values. 

###################################### EXAMPLE USAGE ###################################################

python3 contact_maps.py -i ../../datasets/pdb_structures/ -o ../../datasets/cmaps/ -c 10 -d 10
python3 contact_maps.py -i ../../datasets_frag/pdb_structures/ -o ../../datasets_frag/cmaps_12a/ -c 10 -d 12

WARNING: requires good amount of memory

########################################################################################################
"""


import pandas as pd
import os
import numpy as np
import argparse
import multiprocessing
import matplotlib.pyplot as plt
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB.PDBParser import PDBParser
import warnings
from utils import get_pdb_path_list_from_dir, get_af_id_from_pdb



def get_structure(file):
    """
    Use Bio.PDB to parse protein structure.
    """
    id = file.split("AF-")[1].split("-model_v1")[0]
    return PDBParser().get_structure(id=id, file = f"{file}")[0]


def get_dist_matrix(chain) :
    """
    Compute the distance matrix between C-alpha of a protein.
    """
    m = np.zeros((len(chain), len(chain)), float)
    for i, res1 in enumerate(chain) :
        for j, res2 in enumerate(chain) :
            m[i, j] = abs(res1["CA"] - res2["CA"])
    return m


def get_contact_map(chain, distance=10):
    """
    Compute the contact map between C-alpha of a protein.
    """
    dist_matrix = get_dist_matrix(chain)
    return dist_matrix < distance


def get_3to1_protein_id(protein_id):
    """
    Convert a 3 letter protein code into 1 letter.
    """
    return protein_letters_3to1[protein_id.lower().capitalize()]


def plot_contact_map(contact_map, title=False):
    """
    Plot a 2D 20x20 contact map from contact map matrix.
    """
    aa = [a for a in "ACDEFGHIKLMNPQRSTVWY"]
    im = plt.imshow(contact_map)
    cbar = plt.colorbar(im)
    cbar.set_label("Contact within 10Å")
    plt.yticks(range(20), aa, fontsize=9)
    plt.xticks(range(20), aa, fontsize=9)
    plt.xlabel("AA")
    plt.ylabel("AA")
    if title != False:
        plt.title(title)
    plt.show()


def get_contact_maps(files, output_path, distance=10, verbose=False, num_process=0):
    """
    Given a list of path of PDB file, compute the contact map of each PDB 
    structure and save it as individual .npy file in the given output path.
    """

    # Iterate through the files and save contact map
    for n, file in enumerate(files):
        identifier = get_af_id_from_pdb(file)
        try:
            structure = get_structure(file)
            contact_map = get_contact_map(structure["A"], distance=distance)        
            np.save(f"{output_path}/{identifier}.npy", contact_map.astype(int))
        except:
            warnings.warn(f"WARNING........... could not process {identifier}")
            with open(f"{output_path}/ids_not_processed.txt", 'a+') as file:
                file.write(identifier + '\n')

        # Monitor processing
        if verbose and n % 100 == 0:
            if n == 0:
                print(f"Process [{num_process}] starting..")
            elif n == len(files):
                print(f"Process [{num_process}] completed [{n}/{len(files)}] structures")
            else:
                print(f"Process [{num_process}] completed [{n}/{len(files)}] structures")


def main():

    ## Parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input directory with PDB structures", type=str, required=True)
    parser.add_argument("-o", "--output", help="Output directory to save contact maps", type=str, default="../../datasets/cmaps/")
    parser.add_argument("-d", "--distance", help="Set the distance in angstrom to define a contact", type=int, default=10)
    parser.add_argument("-c", "--num_cores", help="Set the number of cores for parallel processing", type=int)
    parser.add_argument("-v", "--verbose", help="Verbose", type=int, default=1)

    args = parser.parse_args()
    distance = args.distance
    input = args.input
    output = args.output
    num_cores = args.num_cores
    if num_cores is None:
        num_cores = multiprocessing.cpu_count()
    verbose = args.verbose
    print("\nComputing contact maps of all structures in directory")
    print("\nInput directory:", input)
    print("Output:", output)
    print("Distance:", distance)
    print("Num cores:", num_cores)
    print("Verbose:", verbose, "\n")

    # Create necessary folder
    if not os.path.exists(output):
        os.makedirs(output)

    # Get the path of all pdb files in the directorys
    pdb_path_lst = get_pdb_path_list_from_dir(input)

    # Split the PDB files into chunks for each process
    chunk_size = int(len(pdb_path_lst) / num_cores) + 1
    chunks = [pdb_path_lst[i : i + chunk_size] for i in range(0, len(pdb_path_lst), chunk_size)]

    # Create a pool of processes and compute the cmaps in parallel
    with multiprocessing.Pool(processes = num_cores) as pool:
        results = pool.starmap(get_contact_maps, [(chunk, output, distance, verbose, n) 
                                                  for n, chunk in enumerate(chunks)])

    print(f"\nIndividual cmaps saved in {output}")


if __name__ == "__main__":
    main()
