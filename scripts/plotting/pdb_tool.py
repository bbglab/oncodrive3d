import logging
import os

import daiquiri
import subprocess
import pandas as pd
import numpy as np
import re
from tqdm import tqdm

from scripts import __logger_name__

logger = daiquiri.getLogger(__logger_name__ + ".annotations.pdb_tool")

logging.getLogger('urllib3.connectionpool').setLevel(logging.WARNING)



#=========
# pdb_tool
#=========

def run_pdb_tool(pdb_tool_sif_path, i, o, f="4"):
    """
    Use PDB_Tool to extract features from all pdb files in directory.
    """
    
    pdb_files = [file for file in os.listdir(i) if file.endswith(".pdb")]    
    for file in tqdm(pdb_files, desc="Running PDB_Tool"):
        output = f"{o}/{file}".replace(".pdb", ".feature")
        run = subprocess.run(["singularity", "exec", f"{pdb_tool_sif_path}", "/PDB_Tool/PDB_Tool", "-i", f"{i}/{file}", "-o", output, "-F", f])
            
            
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
        df = load_pdb_tool_file(f"{path}/{file}")
        identifier = file.split("-")
        print(identifier)
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