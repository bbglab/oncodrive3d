import os 
import re
import pandas as pd
import subprocess
import argparse
from progressbar import progressbar
import shutil
import warnings
from Bio.PDB import PDBParser

"""
In-house script that merge all fragmented structures in a given directory.
To merge the pdb structures it uses the DEGRONOPEDIA script.

Example usage:

python3 af_merge_parser.py -i ./example_uncompresseds
python3 af_merge_parser.py -i /workspace/projects/clustering_3d/clustering_3d/datasets_frag/pdb_structures
python3 af_merge_parser.py -i ./example_uncompressed -o ./example_uncompressed/out/ 

/workspace/projects/clustering_3d/clustering_3d/datasets_frag/pdb_structures
"""


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
        record = 'SEQRES {:>3} {} {:>4}  {:52}\n'.format(record_counter + 1, "A", num_residues, ' '.join(residue_subset))
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
    Given a directory including pdb files, return 
    a list of tuples (Uniprot_ID, max AF_F).
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


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_dir", help = "Path to directory including pdb structures", type=str, required=True)
parser.add_argument("-o", "--output_dir", help = "Path to output directory", type=str)
parser.add_argument("-v", "--af_version", help = "AlphaFold 2 version used to produced the pdb files", type=int, default=4)
parser.add_argument("-g", "--zip", help = "1 if files is compressed with gzip, else 0", type=int, default=0)

args = parser.parse_args()
input_dir = args.input_dir
output_dir = args.output_dir
af_version = args.af_version
zip = args.zip

dir_path = os.path.abspath(os.path.dirname(__file__))
path_script = f"{dir_path}/af_merge.py"

if output_dir is None:
    output_dir = input_dir
    
if zip:
    zip_ext = ".gzip"
else:
    zip_ext = ""
    
# Create dir where to move original fragmented structures    
path_original_frag = f"{output_dir}/fragmented_pdbs/"
if not os.path.exists(path_original_frag):
    os.makedirs(path_original_frag)
    
# Get list of fragmented Uniprot ID and max AF-F
not_processed = []
fragments = get_list_fragmented_pdb(input_dir)
print("Merging structures..")
for uni_id, max_f in progressbar(fragments):
    
    merge_command = [f"python3", path_script, "-u", uni_id, 
                     "-i", input_dir, "-v", f"{af_version}", "-o", output_dir, "-g", f"{zip}"]
    processed = 0
    
    try:
        run = subprocess.run(merge_command, check=True)
        processed = 1
    except:
        warnings.warn(f"WARNING........... could not process {uni_id} ({max_f} fragments)")
        not_processed.append(uni_id)
        os.remove(f"{output_dir}/AF-{uni_id}-FM-model_v{af_version}.pdb")
    
    # Move the original fragmented structures
    for f in range(1, max_f+1):
        file = f"AF-{uni_id}-F{f}-model_v{af_version}.pdb{zip_ext}"
        shutil.move(f"{input_dir}/{file}", path_original_frag)
        
    # Rename merged structure and add refseq records to pdb
    if processed:
        tmp_name = f"{output_dir}/AF-{uni_id}-FM-model_v{af_version}.pdb"
        name = f"{output_dir}/AF-{uni_id}-F{max_f}M-model_v{af_version}.pdb"
        os.rename(tmp_name, name)
        add_refseq_record_to_pdb(name)

save_unprocessed_ids(not_processed, f"{output_dir}/fragmented_pdbs/not_processed.txt")
