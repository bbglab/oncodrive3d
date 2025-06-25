import os
import gzip
import shutil
import daiquiri

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


def parse_fasta_to_three_letter(path_fasta):
    """
    Parse a FASTA file and return a list of three-letter residue codes.

    :param path_fasta: Path to FASTA file containing single-letter sequence
    :return: List of three-letter codes or raises KeyError on unknown letter
    """
    seq = ''
    with open(path_fasta) as fh:
        for line in fh:
            if line.startswith('>'): continue
            seq += line.strip()
    return [one_to_three_res_map[aa] for aa in seq]


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
    fasta_dir: str = None
    ) -> None:
    """
    Copy all .pdb files from src_dir to dst_dir, renaming and gzipping them.
    Optionally add residues information (SEQRES) from fasta file. 
    """
    
    # Ensure destination directory exists
    os.makedirs(dst_dir, exist_ok=True)

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
        if fasta_dir is not None:
            fasta_path = os.path.join(fasta_dir, f"{accession}.fasta")
            seq = parse_fasta_to_three_letter(fasta_path)
            add_seqres_to_pdb(dst_path, seq)
            logger.debug(f'Inserted SEQRES records into: {new_name}')
        else:
            logger.warning(f'FASTA not found for custom accession {accession}: {fasta_path}')
        
        