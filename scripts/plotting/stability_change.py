import os
import json
import logging
import re
import glob
from multiprocessing import Pool

import pandas as pd
import numpy as np
from progressbar import progressbar
import daiquiri

from scripts import __logger_name__
from scripts.datasets.utils import download_single_file, extract_zip_file
from scripts.globals import rm_dir

logger = daiquiri.getLogger(__logger_name__ + ".plotting.stability_change")

logging.getLogger('urllib3.connectionpool').setLevel(logging.WARNING)


# UniProt accession regex covering both the 6-char (P12345) and
# extended 10-char (A0A075B6X5) forms. See
# https://www.uniprot.org/help/accession_numbers
_UNIPROT_RE = re.compile(
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]"
    r"|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}"
)

# Matches a single missense variant (e.g. "M1A") in RaSP-style CSV rows.
_VARIANT_RE = re.compile(r"([A-Za-z])(\d+)([A-Za-z])")


# ===============================
# Stability change upon mutations
# ===============================


def download_stability_change(path: str,
                              threads: int = 1):
    """
    Downloads stability change upon mutations predicted on AlphaFold 
    structures by RaSP.
    
    Rapid protein stability prediction using deep learning representations
    https://elifesciences.org/articles/82593
    DOI: 10.7554/eLife.82593
    """

    url_website = "https://sid.erda.dk/cgi-sid/ls.py?share_id=fFPJWflLeE"
    filename = "rasp_preds_alphafold_UP000005640_9606_HUMAN_v2.zip"
    download_url = "https://sid.erda.dk/share_redirect/fFPJWflLeE/rasp_preds_alphafold_UP000005640_9606_HUMAN_v2.zip"

    logger.debug(f"Filename: {filename}")
    logger.debug(f"Website url: {url_website}")
    file_path = os.path.join(path, filename)

    try:
        # Download file
        logger.debug(f'Downloading to {file_path}')
        download_single_file(download_url, file_path, min(threads, 4))
        
        # Extract from zip
        logger.debug(f'Extracting {filename}')
        extract_zip_file(file_path, path)
        if os.path.exists(file_path): 
            logger.debug(f'rm {file_path}')
            os.remove(file_path)                       

        logger.debug('Download stability change: SUCCESS')
        logger.debug(f"Files downloaded in directory {path}")
        
        return file_path.replace(".zip", "")

    except Exception as e:
        logger.error('Download stability change: FAIL')
        logger.error(f"Error while downloading stability change: {e}")
        raise e


def append_ddg_to_dict(ddg_dict, df, frag=False):

    pattern = re.compile(r'([A-Za-z])(\d+)([A-Za-z])')
    
    for _, row in df.iterrows():
        variant, ddg = row.values
        pos, alt = extract_mut(variant, pattern)
        
        if pos not in ddg_dict:
            ddg_dict[pos] = {}
        
        if alt not in ddg_dict[pos] and frag:
            ddg_dict[pos][alt] = []

        if frag:
            ddg_dict[pos][alt].append(ddg)
        else:
            ddg_dict[pos][alt] = ddg

    return ddg_dict


def extract_mut(variant_str, pattern):

    match = pattern.match(variant_str)
    pos = match.group(2)
    alt = match.group(3)

    return pos, alt


def save_json(path_dir, uni_id, dictionary):
    
    with open(os.path.join(path_dir, f"{uni_id}_ddg.json"), "w") as json_file:
        json.dump(dictionary, json_file)


def id_from_ddg_path(path):
    """Extract a UniProt accession from a ΔΔG filename.

    Matches the canonical UniProt accession pattern anywhere in the filename,
    so any separator (``-``, ``_``, etc.) around the accession is tolerated.
    Raises ``ValueError`` if no accession is found.
    """

    match = _UNIPROT_RE.search(os.path.basename(path))
    if match is None:
        raise ValueError(
            f"Could not extract a UniProt accession from filename: {path}"
        )
    return match.group(0)


def _validate_protein_ddg(canonical_seq, csv_paths, wt_mismatch_threshold):
    """Validate that predicted ΔΔG variants align with the canonical sequence.

    Two failure modes are caught:
      1. positions outside the canonical sequence (length mismatch / wrong isoform);
      2. wild-type residues disagreeing with the canonical sequence above
         ``wt_mismatch_threshold`` (fraction of total variants checked).

    ``canonical_seq`` is the canonical UniProt sequence string for this protein.
    When ``None`` validation is treated as a hard failure (the caller has
    already decided this protein has no canonical reference available).

    Returns ``(is_valid, reason)`` where ``reason`` is a short string for logging.
    """

    if canonical_seq is None:
        return False, "uniprot_id_not_in_canonical_seq_map"
    if not isinstance(canonical_seq, str) or len(canonical_seq) == 0:
        return False, "canonical_seq_unavailable"

    n_total = 0
    n_mismatch = 0
    for path_prot in csv_paths:
        try:
            df = pd.read_csv(path_prot, usecols=["variant"])
        except Exception as e:
            return False, f"csv_unreadable: {type(e).__name__}"
        for variant in df["variant"]:
            m = _VARIANT_RE.match(str(variant))
            if not m:
                continue
            wt = m.group(1).upper()
            pos = int(m.group(2))
            if pos < 1 or pos > len(canonical_seq):
                return False, (
                    f"position_{pos}_out_of_range_(canonical_len_{len(canonical_seq)})"
                )
            n_total += 1
            if canonical_seq[pos - 1].upper() != wt:
                n_mismatch += 1

    if n_total == 0:
        return False, "no_valid_variants_in_csv"

    mismatch_rate = n_mismatch / n_total
    if mismatch_rate > wt_mismatch_threshold:
        return False, (
            f"wt_mismatch_rate_{mismatch_rate:.1%}_exceeds_threshold_"
            f"{wt_mismatch_threshold:.1%}"
        )
    return True, "ok"


def parse_ddg_rasp_worker(args):

    file, path_dir, output_path, canonical_seq, wt_mismatch_threshold, validate = args

    # Get Uniprot_ID
    try:
        uni_id = id_from_ddg_path(file)
    except ValueError as e:
        logger.warning(f"Skipping ΔΔG file {file}: {e}")
        return

    # Get paths of all fragments
    lst_path_prot = glob.glob(os.path.join(path_dir, f"*{uni_id}*"))
    frag = True if len(lst_path_prot) > 1 else False

    # Validate against the canonical sequence (skipped only if explicitly disabled)
    if validate:
        is_valid, reason = _validate_protein_ddg(
            canonical_seq, lst_path_prot, wt_mismatch_threshold
        )
        if not is_valid:
            logger.warning(f"Skipping ΔΔG for {uni_id}: {reason}")
            return

    # Save a dictionary for each pos with keys as ALT and lst of DDG as values
    ddg_dict = {}
    for path_prot in progressbar(lst_path_prot):
        df = pd.read_csv(path_prot)[["variant", "score_ml"]]
        ddg_dict = append_ddg_to_dict(ddg_dict, df, frag=frag)

    # Iterate through the pos and the ALT and get the mean across frags for each variant
    if frag:
        for pos in ddg_dict:
            for alt in ddg_dict[pos]:
                ddg_dict[pos][alt] = np.mean(ddg_dict[pos][alt])

    # Save dict
    save_json(output_path, uni_id, ddg_dict)


def parse_ddg_rasp(input_path, output_path, threads=1,
                   seq_map=None, wt_mismatch_threshold=0.1):
    """
    It iterates through the csv files in <path_dir> and convert each one into
    a .json dictionary of dictionaries having protein position as keys (str) and
    ALT amino acid (1-letter) as sub-dictionaries keys whose values are the DDG
    (protein stability change upon mutations) for each variant predicted by RaSP.
    If a the protein is fragmented, the DDG of a variant is computed as average
    DDG of that variant across the different fragments (fragments are overlapping).

    When ``seq_map`` is provided (``{uniprot_id: canonical_sequence}``), each
    protein's predictions are validated against its canonical sequence: any
    protein with positions beyond the sequence length or with a wild-type
    mismatch rate above ``wt_mismatch_threshold`` (default 0.1) is skipped
    with a warning. Setting ``wt_mismatch_threshold`` to ``1.0`` effectively
    disables the WT-mismatch check (the position-out-of-range check remains).
    The plotting module already tolerates per-protein missing JSONs, so
    downstream output simply omits the ΔΔG track for skipped genes.

    Only the canonical sequence string for the protein being processed is
    passed to each worker — not the full ``seq_map`` — to keep
    multiprocessing pickle overhead bounded.

    Rapid protein stability prediction using deep learning representations
    https://elifesciences.org/articles/82593
    DOI: 10.7554/eLife.82593
    """

    validate = seq_map is not None
    if seq_map is None:
        seq_map = {}

    # Get already processed files and available ones for processing.
    # Collect (file, uni_id) so we don't extract the accession twice.
    files_processed = set(glob.glob(os.path.join(output_path, "*.json")))
    lst_files = []
    for file in os.listdir(input_path):
        if not file.endswith(".csv"):
            continue
        try:
            uni_id = id_from_ddg_path(file)
        except ValueError as e:
            logger.warning(f"Skipping {file}: {e}")
            continue
        if os.path.join(output_path, f"{uni_id}_ddg.json") in files_processed:
            continue
        lst_files.append((file, uni_id))
    ## Save dict for each proteins
    logger.debug(f"Input: {input_path}")
    logger.debug(f"Output: {output_path}")
    if len(lst_files) > 0:
        logger.debug(f"Parsing DDG of {len(lst_files)} proteins...")

        # TODO: for now it is created a process for each protein, while it would
        #       be better to have chunks of protein processed by the same process
        #       to decrese runtime (at the moment quite slow, 1h40m with 40 cores)

        # TODO: also the parsing itself can be optimized

        # Create a pool of workers parsing processes. Pass only this protein's
        # canonical sequence to each task (avoids pickling the full seq_map
        # ~thousands of times).
        with Pool(processes=threads) as pool:
            args_list = [
                (file, input_path, output_path,
                 seq_map.get(uni_id), wt_mismatch_threshold, validate)
                for file, uni_id in lst_files
            ]
            # Map the worker function to the arguments list
            pool.map(parse_ddg_rasp_worker, args_list)
        if len(lst_files) > 50:
            os.system('clear')
            logger.debug("clear")
        logger.debug("DDG succesfully converted into json files...")
    else:
        logger.debug("DDG not found: Skipping...")

    # Remove the original folder
    logger.debug(f"Deleting {input_path}")
    rm_dir(input_path)
    logger.info("Parsing of DDG completed!")