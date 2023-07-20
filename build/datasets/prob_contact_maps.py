"""
Compute contact maps related features from AlphaFold predicted structures.
The output can be either individual contact maps named with UniprotID_F 
(.npy or .csv), where F is the number of AF fragment. It can generate, or
a single dictionary having UniprotID_F as keys and contact maps as values. 

###################################### EXAMPLE USAGE ###################################################

python3 contact_maps.py -i ../../datasets/pdb_structures/ -o ../../datasets/cmaps/ -c 10 -d 10
    python3 prob_contact_maps.py -i ../../datasets_frag/pdb_structures/ -p ../../datasets_frag/pae/ -o ../../datasets_frag/prob_cmaps/ -c 10

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
import re
import warnings
from math import pi
from utils import get_pdb_path_list_from_dir, get_af_id_from_pdb


# Functions to compute the probability of contact

def cap(r, h):
    """
    Volume of the polar cap of the sphere with radius r and cap height h
    """
    return pi * (3 * r - h) * h ** 2 / 3

def vol(r):
    """
    Volume of sphere with radius r
    """
    return 4 * pi * r ** 3 / 3

def s2_minus_s1(r1, r2, d):
    
    """
    Volume of S2 outside of S1
    r1: radius of S1
    r2: radius of S2
    d: distance between center of S1 and center of S2
    """
    
    # S1 and S2 not in contact
    if d > r1 + r2:      
        return vol(r2)    
    
    # S1 sits inside S2
    elif (r2 > d + r1):  
        return vol(r2) - vol(r1)
        
    # S2 sits inside S1 
    elif (r1 > d + r2):  
        return 0.
    
    # Center of S2 is outside S1
    elif d > r1:         
        h1 = (r2 ** 2 - (d - r1) ** 2) / (2 * d)
        h2 = r1 - h1 + r2 - d
        return vol(r2) - cap(r1, h1) - cap(r2, h2)
    
    # Center of S2 is inside S1
    elif d <= r1:         
        h1 = (r2 ** 2 - (r1 - d) ** 2) / (2 * d)
        h2 = d - r1 + h1 + r2
        return cap(r2, h2) - cap(r1, h1)
    
    
# Other functions

def get_structure(file):
    """
    Use Bio.PDB to parse protein structure.
    """
    
    id = file.split("AF-")[1].split("-model_v1")[0]
    
    return PDBParser().get_structure(id=id, file = f"{file}")[0]


def get_3to1_protein_id(protein_id):
    """
    Convert a 3 letter protein code into 1 letter.
    """
    
    return protein_letters_3to1[protein_id.lower().capitalize()]


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
    dist_matrix = dist_matrix < distance
    
    return dist_matrix.astype(float)



def get_prob_contact(pae_value, dmap_value, distance=10):
    """
    Get probability of contact considering the distance 
    between residues in the predicted structure and the Predicted Aligned Error (PAE).
    """
    if pae_value == 0 and dmap_value == 0:
        
        return 1
    
    else:
    
        # Get the volume of res2 outside of res1
        vol_s2_out_s1 = s2_minus_s1(r1=distance, r2=pae_value, d=dmap_value)

        # Get the probability that s2 is out of s1
        p_s2_in_s1 = vol_s2_out_s1 / vol(pae_value)

        return 1 - p_s2_in_s1
    
    
def get_prob_cmap(chain, pae, distance=10) :
    """
    Compute the probabilities that each pair of residue in a protein are 
    in contact taking into account the Predicted Aligned Error (PAE) and 
    the PDB structure predicted by AlphaFold 2
    """
    
    m = np.zeros((len(chain), len(chain)), float)
    
    for i, res1 in enumerate(chain):
        for j, res2 in enumerate(chain):
            d = abs(res1["CA"] - res2["CA"])
            m[i, j] = get_prob_contact(pae[i, j], d, distance)
            
    return m

def get_prob_cmaps(pdb_files, pae_path, output_path, distance=10, verbose=False, num_process=0):
    """
    Given a list of path of PDB file, compute the probabilistic cmap of 
    each PDB non-fragmented structure and save it as individual .npy file 
    in the given output path. For fragmented structures simply get cmaps.
    """

    # Iterate through the files and save probabilsitic cmap
    for n, file in enumerate(pdb_files):
        identifier = get_af_id_from_pdb(file)
        
        # Get fragmented number
        af_f = identifier.split("-F")[1]
        if af_f.isnumeric():
            af_f = int(af_f)
        else:   
            af_f = int(re.sub(r"\D", "", af_f))
                
        # Get CMAP for fragmented structures (PAE not available yet)    
        if af_f > 1:
            try:
                cmap = get_contact_map(get_structure(file)["A"], distance=distance)    
                np.save(f"{output_path}/{identifier}-prob_cmap.npy", cmap)
            except:
                warnings.warn(f"WARNING........... could not process {identifier}")
                with open(f"{output_path}/ids_not_processed.txt", 'a+') as file:
                    file.write(identifier + '\n')
                    
        # Get probabilistic CMAP
        else:   
            try:
                pae = np.load(f"{pae_path}/{identifier}-predicted_aligned_error.npy")
                chain = get_structure(file)["A"]
                prob_cmap = get_prob_cmap(chain, pae, distance=distance)        
                np.save(f"{output_path}/{identifier}-prob_cmap.npy", prob_cmap)
            except:
                warnings.warn(f"WARNING........... could not process {identifier}")
                with open(f"{output_path}/ids_not_processed.txt", 'a+') as file:
                    file.write(identifier + '\n')

        # Monitor processing
        if verbose and n % 100 == 0:
            if n == 0:
                print(f"Process [{num_process}] starting..")
            elif n == len(pdb_files):
                print(f"Process [{num_process}] completed [{n}/{len(pdb_files)}] structures")
            else:
                print(f"Process [{num_process}] completed [{n}/{len(pdb_files)}] structures")


def main():

    ## Parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_pdb", help="Input directory with PDB structures", type=str, required=True)
    parser.add_argument("-p", "--input_pae", help="Input directory with PAE files", type=str, required=True)
    parser.add_argument("-o", "--output", help="Output directory to save prob cmaps", type=str, default="../../datasets/prob_cmaps/")
    parser.add_argument("-d", "--distance", help="Set the distance in angstrom to define a contact", type=int, default=10)
    parser.add_argument("-c", "--num_cores", help="Set the number of cores for parallel processing", type=int)
    parser.add_argument("-v", "--verbose", help="Verbose", type=int, default=1)

    args = parser.parse_args()
    distance = args.distance
    input_pdb = args.input_pdb
    input_pae = args.input_pae
    output = args.output
    num_cores = args.num_cores
    if num_cores is None:
        num_cores = multiprocessing.cpu_count()
    verbose = args.verbose
    print("\nComputing prob cmaps of all structures in directory")
    print("\nInput PDB directory:", input_pdb)
    print("Input PAE directory:", input_pae)
    print("Output:", output)
    print("Distance:", distance)
    print("Num cores:", num_cores)
    print("Verbose:", verbose, "\n")

    # Create necessary folder
    if not os.path.exists(output):
        os.makedirs(output)

    # Get the path of all pdb files in the directorys
    pdb_path_lst = get_pdb_path_list_from_dir(input_pdb)

    # Split the PDB files into chunks for each process
    chunk_size = int(len(pdb_path_lst) / num_cores) + 1
    chunks = [pdb_path_lst[i : i + chunk_size] for i in range(0, len(pdb_path_lst), chunk_size)]

    # Create a pool of processes and compute the cmaps in parallel
    with multiprocessing.Pool(processes = num_cores) as pool:
        results = pool.starmap(get_prob_cmaps, [(chunk, input_pae, output, distance, verbose, n) 
                                                  for n, chunk in enumerate(chunks)])

    print(f"\nIndividual cmaps saved in {output}")


if __name__ == "__main__":
    main()
