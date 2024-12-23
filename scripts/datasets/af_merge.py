"""
Module to merge overlapping fragments produced by AlphaFold 2
for the predictions of proteins larger than 2700 amino acids.

The module uses an adapted version of the code written by the
authors of DEGRONOPEDIA (Natalia A. Szulc, nszulc@iimcb.gov.pl).
DEGRONOPEDIA - a web server for proteome-wide inspection of degrons
doi: 10.1101/2022.05.19.492622.
"""

import gzip
import logging
import os
import re
import shutil
import subprocess
from os import sep

import daiquiri
import pandas as pd
from Bio.PDB import PDBExceptions, PDBParser, Structure, Superimposer
from Bio.PDB.PDBIO import PDBIO
from tqdm import tqdm

from scripts import __logger_name__

logger = daiquiri.getLogger(__logger_name__ + ".build.af_merge")

daiquiri.getLogger('py.warnings').setLevel(logging.ERROR)


## DEGRONOPEDIA script

def degronopedia_af_merge(struct_name, input_path, afold_version, output_path, zip):
    """
    DEGRONOPEDIA script to merge any AlphaFold fragments into a unique structure.

    DEGRONOPEDIA - a web server for proteome-wide inspection of degrons
    doi: 10.1101/2022.05.19.492622
    https://degronopedia.com/degronopedia/about
    """

    ###  ARGUMENTS PARSING  ###

    if input_path[-1] == sep:
        input_path = input_path[:-1]

    if output_path:
        if output_path[-1] == sep:
            output_path = output_path[:-1]
        save_path = output_path
    else:
        save_path = input_path

    ### FIND OUT HOW MANY PIECES ###

    how_many_pieces = 0
    if zip:
        onlyfiles = [f for f in os.listdir(input_path) if f.endswith(".pdb.gz") and os.path.isfile(os.path.join(input_path, f))]
    else:
        onlyfiles = [f for f in os.listdir(input_path) if f.endswith(".pdb") and os.path.isfile(os.path.join(input_path, f))]

    for f in onlyfiles:
        if struct_name in f and f[0] != '.': # do not include hidden files
            how_many_pieces += 1

    ### MERGING ###

    Bio_parser = PDBParser()
    c = 1
    while c < how_many_pieces :

        struct_save_path = os.path.join(save_path, f"AF-{struct_name}-FM-model_v{afold_version}.pdb")

        # Read reference structure
        if c == 1:
            struct_ref_path = os.path.join(input_path, f"AF-{struct_name}-F{c}-model_v{afold_version}.pdb")
            if zip:
                with gzip.open(f'{struct_ref_path}.gz', 'rt') as handle:
                    struct_ref = Bio_parser.get_structure("ref", handle)
            else:
                with open(struct_ref_path, 'r') as handle:
                    struct_ref = Bio_parser.get_structure("ref", handle)
        else:
            struct_ref = Bio_parser.get_structure("ref", struct_save_path)
        model_ref = struct_ref[0]

        # Read structure to superimpose
        struct_si_path = os.path.join(input_path, f"AF-{struct_name}-F{c+1}-model_v{afold_version}.pdb")
        if zip:
            with gzip.open(f'{struct_si_path}.gz', 'rt') as handle:
                structure_to_superpose = Bio_parser.get_structure("ref", handle)
        else:
            with open(struct_si_path, 'r') as handle:
                structure_to_superpose = Bio_parser.get_structure("ref", handle)
        model_to_super = structure_to_superpose[0]

        # Append atoms from the nine last residues except for the very last one (it is C-end, has one more atom more)
        model_ref_atoms = []
        for j in range(len(model_ref['A'])-9, len(model_ref['A'])):
            for atom in model_ref['A'][j]:
                model_ref_atoms.append(atom)

        # Append atoms from the 1191-1999 residues which correspond the abovementioned residues from the reference
        model_to_superpose_atoms = []
        for j in range(1191, 1200):
            for atom in model_to_super['A'][j]:
                model_to_superpose_atoms.append(atom)

        # Superimpose
        sup = Superimposer()
        sup.set_atoms(model_ref_atoms, model_to_superpose_atoms)

        # Update coords of the residues from the structure to be superimposed
        sup.apply(model_to_super.get_atoms())

        # Delete last residue (C-end residue, with one atom more) from the reference structure
        model_ref['A'].detach_child((' ', len(model_ref['A']), ' '))

        # Delete first 1199 residues from the superimposed structure
        for i in range(1, 1200):
            model_to_super['A'].detach_child((' ', i, ' '))

        # Renumber residues in the superimposed structure
        # Do it twice as you cannot assign a number to a residue that another residue already has
        tmp_resnums = [i+1 for i in range(len(model_to_super['A']))]

        for i, residue in enumerate(model_to_super['A'].get_residues()):
            res_id = list(residue.id)
            res_id[1] = tmp_resnums[i]
            residue.id = tuple(res_id)

        new_resnums = [i+len(model_ref['A'])+1 for i in range(len(model_to_super['A']))]

        for i, residue in enumerate(model_to_super['A'].get_residues()):
            res_id = list(residue.id)
            res_id[1] = new_resnums[i]
            residue.id = tuple(res_id)

        # Merge and save both structures however as two models
        merged = Structure.Structure("master")
        merged.add(model_ref)
        model_to_super.id='B'
        merged.add(model_to_super)

        io = PDBIO()
        io.set_structure(merged)
        io.save(struct_save_path)

        # Unify models
        bashCommand1 = os.path.join("sed '", "TER", f"d' {save_path}", f"AF-{struct_name}-FM-model_v{afold_version}.pdb > {save_path}", "tmp.pdb")
        bashCommand2 = os.path.join("sed '", "MODEL", f"d' {save_path}", f"tmp.pdb > {save_path}", "tmp1.pdb")
        bashCommand3 = os.path.join("sed '", "ENDMDL", f"d' {save_path}", f"tmp1.pdb > {save_path}", "tmp2.pdb")

        subprocess.run(bashCommand1, check=True, text=True, shell=True)
        subprocess.run(bashCommand2, check=True, text=True, shell=True)
        subprocess.run(bashCommand3, check=True, text=True, shell=True)

        # Re-read the structure in Biopython and save
        structure_ok = Bio_parser.get_structure("ok", os.path.join(save_path, "tmp2.pdb"))
        io.set_structure(structure_ok)
        io.save(struct_save_path)

        c += 1

    # Add MODEL 1 at the beggining of the file
    bashCommand4 = os.path.join(f"sed -i '1iMODEL        1                                                                  ' {save_path}", f"AF-{struct_name}-FM-model_v{afold_version}.pdb")
    subprocess.run(bashCommand4, check=True, shell=True)

    # Provide proper file ending
    bashCommand5 = os.path.join(f"sed -i '$ d' {save_path}", f"AF-{struct_name}-FM-model_v{afold_version}.pdb")
    bashCommand6 = os.path.join(f"echo 'ENDMDL                                                                          ' >> {save_path}", f"AF-{struct_name}-FM-model_v{afold_version}.pdb")
    bashCommand7 = os.path.join(f"echo 'END                                                                             ' >> {save_path}", f"AF-{struct_name}-FM-model_v{afold_version}.pdb")
    subprocess.run(bashCommand5, check=True, shell=True)
    subprocess.run(bashCommand6, check=True, text=True, shell=True)
    subprocess.run(bashCommand7, check=True, text=True, shell=True)

    # Delete tmp files
    bashCommand8 = os.path.join(f"rm {save_path}", f"tmp.pdb {save_path}", f"tmp1.pdb {save_path}", "tmp2.pdb")
    subprocess.run(bashCommand8, check=True, text=True, shell=True)


## In-house scripts

# Add SEQREF record to pdb file

def get_res_from_chain(pdb_path):
    """
    Get sequense of amino acid residues from the structure chain.
    """

    # Load structure
    parser = PDBParser()
    structure = parser.get_structure("ID", pdb_path)

    # Get seq from chain
    residues = []
    chain = structure[0]["A"]
    for residue in chain.get_residues():
        residues.append(residue.resname)

    return residues


def get_pdb_seqres_records(lst_res):
    """
    Construct the fixed-width records of a pdb file.
    """

    records = []
    num_residues = len(lst_res)
    record_counter = 0
    while record_counter * 13 < num_residues:
        start_idx = record_counter * 13
        end_idx = min(start_idx + 13, num_residues)
        residue_subset = lst_res[start_idx:end_idx]
        record = 'SEQRES {:>3} {} {:>4}  {:52}\n'.format(record_counter+1, "A", num_residues, ' '.join(residue_subset))
        records.append(record)
        record_counter += 1

    return records


def add_refseq_record_to_pdb(path_structure):
    """
    Add the SEQREF records to the pdb file.
    """

    # Open the PDB file and get SEQRES insert index
    with open(path_structure, 'r') as file:
        pdb_lines = file.readlines()
        insert_index = next(i for i, line in enumerate(pdb_lines) if line.startswith('MODEL'))

    # Get seares records
    residues = get_res_from_chain(path_structure)
    seqres_records = get_pdb_seqres_records(residues)

    # Insert the SEQRES records
    pdb_lines[insert_index:insert_index] = seqres_records

    # Save
    with open(path_structure, 'w') as output_file:
        output_file.truncate()
        output_file.writelines(pdb_lines)


# Other functions

def get_list_fragmented_pdb(pdb_dir):
    """
    Given a directory including pdb files, return a list of tuples (Uniprot_ID, max AF_F).
    """

    # List pdb files
    list_pdb = os.listdir(pdb_dir)
    list_pdb = [file for file in list_pdb if not file.startswith("tmp") and file.endswith(".pdb") or file.endswith(".pdb.gz")]
    list_pdb = [(file.split("-")[1], re.sub(r"\D", "", file.split("-")[2])) for file in list_pdb if file.split("-")[2][-1] != "M"]

    # Get df with max fragment
    df = pd.DataFrame(list_pdb, columns=["Uniprot_ID", "F"])
    df["F"] = pd.to_numeric(df["F"])
    df = df.groupby("Uniprot_ID").max()

    # Get fragmented structures as list of (Uniprot_ID AF_F) tuples
    df = df[df["F"] > 1].reset_index()

    return list(df.to_records(index=False))


def save_unprocessed_ids(uni_ids, filename):

    with open(filename, 'a') as file:
        for id in uni_ids:
            file.write(id + '\n')


# Wrapper function

def merge_af_fragments(input_dir, output_dir=None, af_version=4, gzip=False):
    """
    Run and parse DEGRONOPEDIA script to merge any AlphaFold fragments into a unique structure.

    DEGRONOPEDIA - a web server for proteome-wide inspection of degrons
    doi: 10.1101/2022.05.19.492622
    https://degronopedia.com/degronopedia/about
    """

    if output_dir is None:
        output_dir = input_dir
    if gzip:
        zip_ext = ".gz"
    else:
        zip_ext = ""

    fragments = get_list_fragmented_pdb(input_dir)
    if len(fragments) > 0:

        # Create dir where to move original fragmented structures
        path_original_frag = os.path.join(output_dir, "fragmented_pdbs")
        if not os.path.exists(path_original_frag):
            os.makedirs(path_original_frag)
        checkpoint = os.path.join(path_original_frag, '.checkpoint.merge.txt')

        if os.path.exists(checkpoint):
            logger.debug("Merge fragments already performed: Skipping...")
        else:
            # Get list of fragmented Uniprot ID and max AF-F
            not_processed = []
            for uni_id, max_f in tqdm(fragments, total=len(fragments), desc="Merging AF fragments"):

                processed = False

                try:
                    degronopedia_af_merge(uni_id, input_dir, af_version, output_dir, gzip)
                    processed = True
                except (PDBExceptions.PDBIOException, PDBExceptions.PDBConstructionException):
                    logger.warning(f"Could not process {uni_id} ({max_f} fragments)")
                    not_processed.append(uni_id)
                    f_path = os.path.join(output_dir, f"AF-{uni_id}-FM-model_v{af_version}.pdb")
                    if os.path.exists(f_path):
                        os.remove(f_path)
                # Move the original fragmented structures
                for f in range(1, max_f+1):
                    file = f"AF-{uni_id}-F{f}-model_v{af_version}.pdb{zip_ext}"
                    shutil.move(os.path.join(input_dir, file), path_original_frag)

                # Rename merged structure and add refseq records to pdb
                if processed:
                    tmp_name = os.path.join(output_dir, f"AF-{uni_id}-FM-model_v{af_version}.pdb")
                    name = os.path.join(output_dir, f"AF-{uni_id}-F{max_f}M-model_v{af_version}.pdb")
                    os.rename(tmp_name, name)
                    add_refseq_record_to_pdb(name)

            if len(not_processed) > 0:
                logger.warning(f"Not processed: {not_processed}")
            with open(checkpoint, "w") as f:
                f.write('')

            save_unprocessed_ids(not_processed,
                                os.path.join(output_dir, "fragmented_pdbs", "ids_not_merged.txt"))
            logger.info("Merge of structures completed!")

    else:
        logger.debug("Nothing to merge: Skipping...")
