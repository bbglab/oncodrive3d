"""
Module to compute maps of probabilities of contact (pCMAPs) from
AlphaFold predicted structures and predicted aligned error (PAE).

The pCMAPs is a dataframe including the probability of contact for
each pair of residue of a given protein. Given a threshold (10Å as
default) to define the contact, the probability that two residues i
and j are in contact is computed considering the distance between i
and j in the PDB structure and their predicted error in the PAE.
"""


import gzip
import daiquiri
import multiprocessing
import os
import re
from math import pi

import numpy as np
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB.PDBParser import PDBParser

from scripts import __logger_name__
from scripts.datasets.utils import (get_af_id_from_pdb,
                                    get_pdb_path_list_from_dir)

logger = daiquiri.getLogger(__logger_name__ + ".build.prob_contact_maps")


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

    if file.endswith('.gz'):
        with gzip.open(file, 'rt') as handle:
            return PDBParser().get_structure(id=id, file=handle)[0]
    else:
        with open(file, 'r') as handle:
            return PDBParser().get_structure(id=id, file=handle)[0]


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

    return dist_matrix < distance


def get_prob_contact(pae_value, dmap_value, distance=10):
    """
    Get probability of contact considering the distance
    between residues in the predicted structure and the
    Predicted Aligned Error (PAE).
    """

    if pae_value == 0:

        if dmap_value < distance:
            return 1
        else:
            return 0

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


def get_prob_cmaps(pdb_files, pae_path, output_path, distance=10, num_process=0):
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
                np.save(os.path.join(output_path, f"{identifier}.npy"), cmap)
            except Exception as e:
                logger.warning(f"Could not process {identifier}")
                logger.warning(f"Error: {e}")
                with open(os.path.join(output_path, "ids_not_processed.txt"), 'a+') as file:
                    file.write(identifier + '\n')

        # Get probabilistic CMAP
        else:
            try:
                pae = np.load(os.path.join(pae_path, f"{identifier}-predicted_aligned_error.npy"))
                chain = get_structure(file)["A"]
                prob_cmap = get_prob_cmap(chain, pae, distance=distance)
                np.save(os.path.join(output_path, f"{identifier}.npy"), prob_cmap)
            except Exception as e:
                logger.warning(f"Could not process {identifier}")
                logger.warning(f"Error: {e}")
                with open(os.path.join(output_path, "ids_not_processed.txt"), 'a+') as file:
                    file.write(identifier + '\n')

        # Monitor processing
        if n % 100 == 0:
            if n == 0:
                logger.debug(f"Process [{num_process}] starting processing [{len(pdb_files)}] structures...")
            else:
                logger.debug(f"Process [{num_process}] completed [{n}/{len(pdb_files)}] structures...")
        elif n+1 == len(pdb_files):
            logger.debug(f"Process [{num_process}] completed!")



def get_prob_cmaps_mp(input_pdb,
                      input_pae,
                      output,
                      distance = 10,
                      num_cores = 1):
    """
    Given a list of path of PDB file, use multiprocessing to compute pCMAPs
    (maps or probabilities of contact between each residues) for each PDB
    non-fragmented structure and save it as individual .npy file in the given
    output path. For fragmented structures simply get cmaps.
    """

    n_structures = len([pdb for pdb in os.listdir(input_pdb) if ".pdb" in pdb])
    logger.debug(f"Input PDB directory: {input_pdb}")
    logger.debug(f"Input PAE directory: {input_pae}")
    logger.debug(f"Output: {output}")
    logger.debug(f"Distance: {distance}")
    logger.debug(f"Cores: {num_cores}")

    # Create necessary folder
    if not os.path.exists(output):
        os.makedirs(output)

    # Get the path of all pdb files in the directorys
    pdb_path_lst = get_pdb_path_list_from_dir(input_pdb)

    # Split the PDB files into chunks for each process
    chunk_size = int(len(pdb_path_lst) / num_cores) + 1
    chunks = [pdb_path_lst[i : i + chunk_size] for i in range(0, len(pdb_path_lst), chunk_size)]

    # Create a pool of processes and compute the cmaps in parallel
    logger.debug(f"Processing [{n_structures}] structures by [{len(chunks)}] processes...")
    with multiprocessing.Pool(processes = num_cores) as pool:
        results = pool.starmap(get_prob_cmaps,
                               [(chunk, input_pae, output, distance, n) for n, chunk in enumerate(chunks)])
