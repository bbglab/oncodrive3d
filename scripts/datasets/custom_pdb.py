import os
import gzip
import shutil
import daiquiri
import pandas as pd
from scripts.datasets.utils import get_seq_from_pdb

from scripts import __logger_name__

logger = daiquiri.getLogger(__logger_name__ + ".build.custom_PDB")


# Mapping single-letter to three-letter codes
one_to_three_res_map = {
    'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY',
    'H':'HIS','I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN',
    'P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL',
    'W':'TRP','Y':'TYR'
}


def get_pdb_seqres_records(lst_res):
    """
    Given a list of three-letter residue names, construct SEQRES records.

    :param lst_res: List of residue names (e.g. ['MET','ALA',...])
    :return: List of formatted SEQRES lines ending with '\n'
    """
    records = []
    num_res = len(lst_res)
    for i in range(0, num_res, 13):
        block = lst_res[i:i+13]
        rec_num = len(records) + 1
        # Atom name block: up to 13 residues separated by spaces, padded to 51 chars
        residues_str = ' '.join(block)
        records.append(
            f"SEQRES {rec_num:>3} A {num_res:>4}  {residues_str:<51}\n"
        )
    return records


def add_seqres_to_pdb(path_pdb: str, residues: list) -> None:
    """
    Insert SEQRES records at the very top of a PDB file (supports gzipped and plain).

    :param path_pdb: Path to the PDB file to modify in-place (can be .pdb or .pdb.gz)
    :param residues: List of three-letter residue codes for the sequence
    """
    # Choose appropriate open functions
    open_in = gzip.open if path_pdb.endswith('.gz') else open
    mode_in = 'rt' if path_pdb.endswith('.gz') else 'r'
    open_out = gzip.open if path_pdb.endswith('.gz') else open
    mode_out = 'wt' if path_pdb.endswith('.gz') else 'w'

    # Read original lines
    with open_in(path_pdb, mode_in) as fh:
        lines = fh.readlines()

    # Generate SEQRES lines
    seqres = get_pdb_seqres_records(residues)

    # Insert at top
    new_lines = seqres + lines

    # Write back
    with open_out(path_pdb, mode_out) as fh:
        fh.writelines(new_lines)
        
        
def copy_and_parse_custom_pdbs(
    src_dir: str,
    dst_dir: str,
    af_version: int = 4,
    custom_mane_metadata_path: str = None
    ) -> None:
    """
    Copy all .pdb files from src_dir to dst_dir, renaming and gzipping them.
    Optionally add residues information (SEQRES) from samplesheet. 
    """
    
    # Ensure destination directory exists
    os.makedirs(dst_dir, exist_ok=True)
    
    # Handle metadata
    if custom_mane_metadata_path and os.path.isfile(custom_mane_metadata_path):
        samplesheet_df = pd.read_csv(custom_mane_metadata_path)
    else:
        if not custom_mane_metadata_path:
            logger.warning(
                "No samplesheet path provided; skipping SEQRES insertion for accession %s",
                accession
            )
        else:
            logger.warning(
                "samplesheet not found at '%s'; skipping SEQRES insertion for accession %s",
                custom_mane_metadata_path,
                accession
            )
        samplesheet_df = None

    # Copy and gzip pdb and optionally add REFSEQ
    for fname in os.listdir(src_dir):
        if not fname.endswith('.pdb'):
            continue

        parts = fname.split('.')  # e.g. [ACCESSION, fragment_code, 'alphafold', 'pdb']
        if len(parts) < 4:
            logger.warning(f"Skipping unexpected filename format: {fname}")
            continue

        accession = parts[0]
        fragment = parts[1] 

        new_name = f'AF-{accession}-F{fragment}-model_v{af_version}.pdb.gz'

        src_path = os.path.join(src_dir, fname)
        dst_path = os.path.join(dst_dir, new_name)

        with open(src_path, 'rb') as fin, gzip.open(dst_path, 'wb') as fout:
            shutil.copyfileobj(fin, fout)

        logger.debug(f'Copied and gzipped: {fname} -> {new_name}')
        
        # Optionally add SEQRES records
        if samplesheet_df is not None:
            if accession not in samplesheet_df["sequence"].values:
                logger.warning("Accession %s not in samplesheet %s", accession, custom_mane_metadata_path)
                continue
            seq = samplesheet_df[samplesheet_df["sequence"] == accession].refseq.values[0]
            
            if not pd.isna(seq):
                seq = [one_to_three_res_map[aa] for aa in seq]
                add_seqres_to_pdb(path_pdb=dst_path, residues=seq)
                logger.debug(f"Inserted SEQRES records into: {new_name}")
            else:
                try:
                    seq = "".join(list(get_seq_from_pdb(dst_path)))
                    if not pd.isna(seq):
                        logger.debug(f"SEQRES record already present in the structure: {new_name}")
                    else:
                        logger.warning(f"SEQRES record not in samplesheet and could not extract it from the structure: {new_name}")
                except Exception as e:
                    logger.warning(f"SEQRES record not in samplesheet and could not extract it from the structure: {new_name}")
                    logger.warning(f"Exception captured: {e}")