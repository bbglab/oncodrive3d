import os
import logging
import subprocess
import re
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
from tqdm import tqdm
import daiquiri

from scripts import __logger_name__

logger = daiquiri.getLogger(__logger_name__ + ".plotting.pdb_tool")

logging.getLogger('urllib3.connectionpool').setLevel(logging.WARNING)


# ========
# PDB_Tool
# ========


def decompress_pdb_gz(input_dir):
    """
    Extract all .gz PDB files in directory.
    """
    
    files_in_directory = os.listdir(input_dir)
    gz_file_present = any(file.endswith(".pdb.gz") for file in files_in_directory)
    if gz_file_present:
        logger.debug("Decompressing .gz PDB files...")
        command = f'gunzip -q {os.path.join(input_dir, "*.pdb.gz")}'
        subprocess.run(command, shell=True)
        logger.debug("Decompression complete!")
    

def _run_pdb_tool_one(input_path, output_path, f):
    """Invoke the PDB_Tool binary on a single structure. Captures stderr so failures
    surface as warnings instead of garbling the tqdm bar."""
    result = subprocess.run(
        ["PDB_Tool", "-i", input_path, "-o", output_path, "-F", f],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        logger.warning(
            f"PDB_Tool failed on {os.path.basename(input_path)} "
            f"(rc={result.returncode}): {result.stderr.strip()}"
        )


def run_pdb_tool(input_dir, output_dir, f="4", cores=1):
    """
    Use PDB_Tool to extract features from all pdb files in directory.

    Each invocation is independent, so the calls are dispatched across a thread
    pool. PDB_Tool runs as a subprocess, so the GIL doesn't constrain throughput.
    """

    pdb_tool_output = os.path.join(output_dir, "pdb_tool")
    if not os.path.isdir(pdb_tool_output):
        os.makedirs(pdb_tool_output)
        logger.debug(f'mkdir {pdb_tool_output}')

    decompress_pdb_gz(input_dir)
    pdb_files = [file for file in os.listdir(input_dir) if file.endswith(".pdb")]
    if not pdb_files:
        return pdb_tool_output

    max_workers = max(1, min(cores, len(pdb_files)))
    logger.debug(f"Running PDB_Tool with {max_workers} worker(s)...")

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(
                _run_pdb_tool_one,
                os.path.join(input_dir, file),
                os.path.join(pdb_tool_output, file.replace(".pdb", ".feature")),
                f,
            )
            for file in pdb_files
        ]
        for future in tqdm(as_completed(futures), total=len(futures), desc="Running PDB_Tool"):
            future.result()

    return pdb_tool_output
            
            
def load_pdb_tool_file(path):
    """
    Parse .feature file obtained by PDB_Tool.
    """
    
    with open(path, "r") as f:
        lines = f.readlines()
        lst_lines = []
        for l in lines:
            lst_lines.append(l.strip().split())
   # return lst_lines
    df = pd.DataFrame(lst_lines[1:], columns = lst_lines[0]).iloc[:,:9]
    df = df.drop("Missing", axis=1)
    df = df.rename(columns={"Num" : "Pos"})

    for c in df.columns:
        if c == "Pos" or c == "ACC" or c == "CNa" or c == "CNb":
            data_type = int
        else:
            data_type = float
        try:
            df[c] = df[c].astype(data_type)
        except:
            pass
    return df


def get_pdb_tool_file_in_dir(path):
    """
    Get the list of PDB_Tool .feature files from directory.
    """
    
    list_files = os.listdir(path)
    ix = [re.search('.\.feature$', x) is not None for x in list_files]
    
    return list(np.array(list_files)[np.array(ix)])


def load_all_pdb_tool_files(path):
    """
    Get the features of all proteins in the directory.
    """
    
    df_list = []
    feature_file_list = get_pdb_tool_file_in_dir(path)
    
    for file in tqdm(feature_file_list, desc="Parsing PDB_tool output"):
        df = load_pdb_tool_file(os.path.join(path, file))
        identifier = file.split("-")
        df["Uniprot_ID"] = identifier[1]
        df["F"] = identifier[2].replace("F", "")
        df_list.append(df)

    return pd.concat(df_list).reset_index(drop=True)


def pdb_tool_to_3s_sse(df):
    """
    Reduce secondary structure from 8 to 3 states.
    """

    mapper = {"H":"Helix", 
              "G":"Helix", 
              "I":"Helix", 
              "L":"Coil", 
              "T":"Coil", 
              "S":"Coil", 
              "B":"Coil", 
              "E":"Ladder"}
    df["SSE"] = df["SSE"].map(mapper)

    return df


def parse_pdb_tool(input_dir : str, output_dir : str):
    """
    Parse PDB_Tool .feature files inclued in the path into a unique df.
    """
    
    pdb_tool_df = load_all_pdb_tool_files(input_dir)
    pdb_tool_df = pdb_tool_to_3s_sse(pdb_tool_df)
    pdb_tool_df = pdb_tool_df.drop(columns=["CLE", "ACC", "CNa", "CNb"])
    pdb_tool_df.to_csv(os.path.join(output_dir, "pdb_tool_df.tsv"), sep="\t", index=False)
    try:
        logger.debug(f"Deleting {input_dir}")
        shutil.rmtree(input_dir)
    except OSError as e:
        logger.debug(f"Cold not delete {input_dir}\nError: {e}")