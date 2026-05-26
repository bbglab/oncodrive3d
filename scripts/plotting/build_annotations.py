import os
import logging

import daiquiri
import pandas as pd

from scripts import __logger_name__
from scripts.globals import copy_dir
from scripts.plotting.utils import clean_annot_dir, get_species
from scripts.plotting.stability_change import download_stability_change, parse_ddg_rasp
from scripts.plotting.pdb_tool import run_pdb_tool, parse_pdb_tool
from scripts.plotting.pfam import get_pfam
from scripts.plotting.uniprot_feat import get_uniprot_feat

logger = daiquiri.getLogger(__logger_name__ + ".plotting.build")

logging.getLogger('urllib3.connectionpool').setLevel(logging.WARNING)


# =================
# BUILD ANNOTATIONS 
# ================= 

def get_annotations(data_dir,
                    output_dir,
                    ddg_dir,
                    #path_pdb_tool_sif,
                    organism,
                    cores,
                    ddg_mismatch_threshold=0.1):
    """
    Main function to build annotations to generate annotated plots.
    """

    # Empty directory
    clean_annot_dir(output_dir, 'd')
    seq_df = None  # loaded lazily — only the ΔΔG and Pfam/UniProt steps need it

    # # Get ddG
    species = get_species(organism)
    logger.info("Checking stability change inputs..")
    if species == "Mus musculus" and ddg_dir is None:
        logger.warning(
            "Stability change (ΔΔG) for Mus musculus requires precomputed predictions: "
            "the public RaSP bundle is only hosted for Homo sapiens. Skipping ΔΔG step. "
            "Provide --ddg_dir with mouse RaSP predictions to include ΔΔG tracks."
        )
    elif species in ("Homo sapiens", "Mus musculus"):
        ddg_output = os.path.join(output_dir, "stability_change")
        os.makedirs(ddg_output, exist_ok=True)
        if ddg_dir is not None:
            # Copy ddG from path
            temp_ddg_path = os.path.join(output_dir, "stability_change_temp")
            copy_dir(source_dir=ddg_dir, destination_dir=temp_ddg_path)
        else:
            # Download ddG
            temp_ddg_path = download_stability_change(ddg_output, cores)
        logger.info("Completed!")

        ## TODO: Optimize DDG parsing
        ##       - only one protein is allocated to one process every time
        ##       - a list of proteins should be allocated instead

        # Parsing DDG (with per-protein validation against canonical sequence)
        seq_df = pd.read_table(os.path.join(data_dir, "seq_for_mut_prob.tsv"))
        seq_map = dict(zip(seq_df["Uniprot_ID"], seq_df["Seq"]))
        logger.info("Parsing stability change..")
        parse_ddg_rasp(temp_ddg_path, ddg_output, cores,
                       seq_map=seq_map,
                       wt_mismatch_threshold=ddg_mismatch_threshold)
        logger.info("Parsing completed!")
    else:
        logger.warning(f"Currently, stability change annotation is not available for {species} but only for Homo sapiens: Skipping...")
    
    ## TODO: Enable multiprocessing for PDB_Tool
    ## TODO: Evaluate the possibility of installing PDB_Tool instead of using container
    
    # Run PDB_Tool
    logger.info("Extracting pACC and 2° structure..")
    path_pdb_structure = os.path.join(data_dir, "pdb_structures")
    pdb_tool_output = run_pdb_tool(input_dir=path_pdb_structure, output_dir=output_dir)
    logger.info("Extraction completed!")
    
    # Parse PDB_Tool
    logger.info("Parsing pACC and 2° structures..")
    parse_pdb_tool(input_dir=pdb_tool_output, output_dir=output_dir)
    logger.info("Parsing completed!")
    
    # Get Pfam annotations
    logger.info("Downloading and parsing Pfam..")
    if seq_df is None:
        seq_df = pd.read_table(os.path.join(data_dir, "seq_for_mut_prob.tsv"))
    pfam_df = get_pfam(seq_df = seq_df,
                       output_tsv = os.path.join(output_dir, "pfam.tsv"),
                       organism = species)
    logger.info("Completed!")
    
    # Get Uniprot features
    logger.info("Downloading and parsing Features from Uniprot..")
    get_uniprot_feat(seq_df = seq_df, 
                     pfam_df = pfam_df,
                     output_tsv = os.path.join(output_dir, "uniprot_feat.tsv"))      
    logger.info("Completed!")
    

if __name__ == "__main__":
    get_annotations(data_dir = "./datasets",
                    output_dir = "./annotations",
                    organism = "Homo sapiens",
                    cores = 4)