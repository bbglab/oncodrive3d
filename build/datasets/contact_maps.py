"""
Compute contact maps related features from AlphaFold predicted structures.
The output can be either individual contact maps named with UniprotID_F 
(.npy or .csv), where F is the number of AF fragment. It can generate, or
a single dictionary having UniprotID_F as keys and contact maps as values. 

###################################### EXAMPLE USAGE ###################################################

python3 contact_maps.py -i /workspace/datasets/alphafold_features/AF_homo_sapiens_pred/ \
-o /../../datasets/cmaps/

WARNING: requires good amount of memory

########################################################################################################
"""


import pandas as pd
import os
import re
import numpy as np
import argparse
import pickle
import matplotlib.pyplot as plt
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB.PDBParser import PDBParser
from utils import get_pdb_path_list_from_dir, get_seq_from_pdb, get_af_id_from_pdb



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


def get_contact_maps_from_dir(input_path, output_path, distance=10, 
                              output_dict=False, pandas=False, verbose=True):
    """
    Compute the contact maps of all PDB files contained in a given directory.
    Save each contact map as pandas dataframe or numpy array in the given output 
    directory. Save the contact maps as individual files or as values of a 
    dictionary whose keys are Uniprot ID.
    """
    
    # Get the path of all pdb files in the directory
    pdb_path_list = get_pdb_path_list_from_dir(input_path)
    cmaps_dict = {}
    
    # Iterate through the files and save contact map
    for n, file in enumerate(pdb_path_list):
        identifier = get_af_id_from_pdb(file)
        structure = get_structure(file)
        contact_map = get_contact_map(structure["A"], distance=distance)

        # Save as pd dataframes
        if pandas:
            contact_map = pd.DataFrame(contact_map, 
                                       columns = get_seq_from_pdb(file))
            if output_dict:
                cmaps_dict[identifier] = contact_map
            else:
                contact_map.to_csv(f"{output_path}/{identifier}.csv", index=False)
        # Save as np arrays
        else:
            contact_map = contact_map.astype(int)
            if output_dict:
                cmaps_dict[identifier] = contact_map
            else:
                np.save(f"{output_path}/{identifier}.npy", contact_map)

        # Monitor processing
        if verbose and n % 1000 == 0:
            if n == 0:
                print("Processing structures..")
            else:
                print(f"Completed [{n}/{len(pdb_path_list)}] structures")
        
    if output_dict:
        return cmaps_dict
    else:
        return None


def main():

    ## Parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input directory with PDB structures", type=str, required=True)
    parser.add_argument("-o", "--output", help="Output directory to save contact maps", type=str, default="../../datasets/cmaps/")
    parser.add_argument("-d", "--distance", help="Set the distance in angstrom to define a contact", type=int, default=10)
    parser.add_argument("-p", "--pandas", help="1 get pandas dataframe (.csv), 0 get numpy array (.npy) for individual cmaps", type=int, default=0)

    # Optional parameters to get a single dictionary with AF_uniprot ID as keys and cmaps as values
    parser.add_argument("-D", "--dictionary", help="if 1 get dictionary (.pickle), if 0 get individual cmaps file (.npy or .csv)", type=int, default=0)
    parser.add_argument("-f", "--filename", help="Filename for the contact map dictionary", type=str, default="contact_maps_dict_all_prot.pickle")

    args = parser.parse_args()
    distance = args.distance
    filename = args.filename
    dictionary = args.dictionary
    input = args.input
    output = args.output
    pandas = args.pandas
    print("\nComputing contact maps of all structures in directory")
    print("\nInput directory:", input)
    print("Output:", output)
    print("Distance:", distance)
    print("Dictionary:", dictionary, "\n")

    # Create necessary folder
    if not os.path.exists(output):
        os.makedirs(output)

    # Compute the contact maps
    maps_dict = get_contact_maps_from_dir(input, output, distance=distance, output_dict=dictionary, pandas=pandas)

    if dictionary:
        print("Saving dictionary..")
        pickle.dump(maps_dict, open(f"{output}/{filename}", 'wb'))
        print(f"File saved in {output}/{filename}")
    else:
        print(f"Individual cmaps saved in {output}")


if __name__ == "__main__":
    main()
