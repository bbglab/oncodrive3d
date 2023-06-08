#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Script to merge split pdbfiles from AlphaFold.

**Author:**
Natalia A. Szulc, nszulc@iimcb.gov.pl

If you used this script, please check if our paper on DEGRONOPEDIA has been published (see section Citing at https://degronopedia.com/degronopedia/about), if not, please cite the preprint:

Natalia A. Szulc, Filip Stefaniak, Małgorzata Piechota, Andrea Cappannini, Janusz M. Bujnicki, Wojciech Pokrzywa (2022).
DEGRONOPEDIA - a web server for proteome-wide inspection of degrons
doi: 10.1101/2022.05.19.492622

Example: python3 merge_afold_files.py -u A2VEC9 -i example_uncompressed/ -v 4 -o example_uncompressed/out/ -g 0
'''

import argparse
from Bio.PDB import Structure
from Bio.PDB import *
import gzip
from os import listdir, sep
from os.path import isfile, join
import subprocess

if __name__ == "__main__":

    ###  ARGUMENTS PARSING  ###

    parser = argparse.ArgumentParser(description = '''Tool merging split AlphaFold PDB files.''',
                                     add_help = False,
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter,
                                     conflict_handler = 'resolve')

    required_arguments = parser.add_argument_group('Required arguments')
    required_arguments.add_argument('-u', help='UniProt ID', required=True, metavar='UniProtID', default=argparse.SUPPRESS)
    required_arguments.add_argument('-i', help='input path', required=True, metavar='Input path', default=argparse.SUPPRESS)
    required_arguments.add_argument('-v', help='AlphaFold version', required=True, metavar='AlphaFold version', default=argparse.SUPPRESS)

    optional_arguments = parser.add_argument_group('Optional arguments')
    optional_arguments.add_argument('-o', help='output path', metavar='Output path')
    optional_arguments.add_argument('-g', help='gzip file', type=int, default=1)
    optional_arguments.add_argument('-h', action = 'help', help = 'show this help message and exit')
    optional_arguments.add_argument('--help', '-h', action = 'help', help = 'show this help message and exit')

    args = vars(parser.parse_args())

    struct_name = args['u']
    input_path = args['i']
    afold_version = args['v']
    output_path = args['o']
    zip = args['g']

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
    onlyfiles = [f for f in listdir(input_path) if isfile(join(input_path, f))]

    for f in onlyfiles:
        if struct_name in f and f[0] != '.': # do not include hidden files
            how_many_pieces += 1

    ### MERGING ###

    Bio_parser = PDBParser()
    c = 1
    while c < how_many_pieces :

        # Read reference structure
        if c == 1:
            if zip:  
                with gzip.open(f'{input_path}/AF-{struct_name}-F{c}-model_v{afold_version}.pdb.gz', 'rt') as handle:
                    structure_ref = Bio_parser.get_structure("ref", handle)
            else:
                with open(f'{input_path}/AF-{struct_name}-F{c}-model_v{afold_version}.pdb', 'r') as handle:
                    structure_ref = Bio_parser.get_structure("ref", handle)
        else:
            structure_ref = Bio_parser.get_structure("ref", f"{save_path}/AF-{struct_name}-FM-model_v{afold_version}.pdb")
        model_ref = structure_ref[0]

        # Read structure to superimpose
        if zip:  
            with gzip.open(f'{input_path}/AF-{struct_name}-F{c+1}-model_v{afold_version}.pdb.gz', 'rt') as handle:
                structure_to_superpose = Bio_parser.get_structure("ref", handle)
        else:
            with open(f'{input_path}/AF-{struct_name}-F{c+1}-model_v{afold_version}.pdb', 'r') as handle:
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
        io.save(f"{save_path}/AF-{struct_name}-FM-model_v{afold_version}.pdb")

        # Unify models
        bashCommand1 = f"sed '/TER/d' {save_path}/AF-{struct_name}-FM-model_v{afold_version}.pdb > {save_path}/tmp.pdb"
        bashCommand2 = f"sed '/MODEL/d' {save_path}/tmp.pdb > {save_path}/tmp1.pdb"
        bashCommand3 = f"sed '/ENDMDL/d' {save_path}/tmp1.pdb > {save_path}/tmp2.pdb"

        subprocess.run(bashCommand1, check=True, text=True, shell=True)
        subprocess.run(bashCommand2, check=True, text=True, shell=True)         
        subprocess.run(bashCommand3, check=True, text=True, shell=True)

        # Re-read the structure in Biopython and save
        structure_ok = Bio_parser.get_structure("ok", f'{save_path}/tmp2.pdb')
        io.set_structure(structure_ok)
        io.save(f"{save_path}/AF-{struct_name}-FM-model_v{afold_version}.pdb")

        c += 1

    # Add MODEL 1 at the beggining of the file
    #subprocess.run(f"gsed -i -e '1iMODEL        1                                                                  \' {save_path}/AF-{struct_name}-FM-model_v{afold_version}.pdb", check=True, text=True, shell=True)
    subprocess.run(f"sed -i '1iMODEL        1                                                                  ' {save_path}/AF-{struct_name}-FM-model_v{afold_version}.pdb", check=True, shell=True)

    # Provide proper file ending
    #subprocess.run(f"gsed -i '$ d' {save_path}/AF-{struct_name}-FM-model_v{afold_version}.pdb", check=True, text=True, shell=True)
    subprocess.run(f"sed -i '$ d' {save_path}/AF-{struct_name}-FM-model_v{afold_version}.pdb", check=True, shell=True)
    subprocess.run(f"echo 'ENDMDL                                                                          ' >> {save_path}/AF-{struct_name}-FM-model_v{afold_version}.pdb", check=True, text=True, shell=True)
    subprocess.run(f"echo 'END                                                                             ' >> {save_path}/AF-{struct_name}-FM-model_v{afold_version}.pdb", check=True, text=True, shell=True)


    # Delete tmp files
    subprocess.run(f"rm {save_path}/tmp.pdb {save_path}/tmp1.pdb {save_path}/tmp2.pdb", check=True, text=True, shell=True)
