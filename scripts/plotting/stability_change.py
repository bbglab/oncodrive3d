import os
import json
import logging
import re
import glob
from collections import Counter
from multiprocessing import Pool

import pandas as pd
import numpy as np
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

# Missense variant like "M1A". The ^ anchor is required because str.extract uses re.search.
_VARIANT_RE = re.compile(r"^([A-Za-z])(\d+)([A-Za-z])")

# Skip-reason keys returned by parse_ddg_rasp_worker and the human-readable
# labels used in the INFO summary at the end of parse_ddg_rasp.
_SKIP_REASONS = {
    "no_uniprot_match": "UniProt ID not found in datasets/seq_for_mut_prob.tsv",
    "no_seq_available": "Sequence missing in datasets/seq_for_mut_prob.tsv",
    "out_of_range":     "Position out of range",
    "wt_mismatch":      "WT mismatch above threshold",
    "no_variants":      "No valid variants in CSV",
    "unreadable":       "Unreadable CSV",
    "bad_filename":     "UniProt accession not recognised in filename",
}


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


def _parse_ddg_csv(path_prot):
    """Read one RaSP-style CSV and extract (wt, pos, alt, ddg) arrays.

    Returns ``(wt_arr, pos_arr, alt_arr, ddg_arr)`` for variants whose ``variant``
    field matches ``[A-Za-z]\\d+[A-Za-z]``. Malformed rows are silently dropped,
    consistent with the previous row-by-row parser.

    On read error (missing column, broken CSV, etc.) raises the exception so the
    caller can log a clean ``csv_unreadable`` reason.
    """

    # Read only the two columns we need — RaSP files have ~13 columns by default.
    df = pd.read_csv(path_prot, usecols=["variant", "score_ml"])
    # Vectorised regex extract; rows that don't match yield NaN.
    extracted = df["variant"].astype(str).str.extract(_VARIANT_RE, expand=True)
    valid = extracted.notna().all(axis=1)
    if not valid.any():
        empty = np.array([], dtype=object)
        return empty, empty, empty, np.array([], dtype=float)
    wt_arr = extracted.loc[valid, 0].to_numpy()
    pos_arr = extracted.loc[valid, 1].to_numpy()
    alt_arr = extracted.loc[valid, 2].to_numpy()
    ddg_arr = df.loc[valid, "score_ml"].to_numpy(dtype=float)  # int score_ml would yield np.int64 (not JSON-serializable)
    return wt_arr, pos_arr, alt_arr, ddg_arr


def parse_ddg_rasp_worker(args):
    """Process one protein: read its fragment CSV(s), validate, accumulate, write JSON.

    Each fragment CSV is read exactly once. Validation (position-in-range +
    cumulative WT-mismatch rate) is interleaved with accumulation rather than
    run as a separate pre-pass, eliminating the double-read while preserving
    the original drop-protein-on-failure semantics.
    """

    file, path_dir, output_path, canonical_seq, wt_mismatch_threshold, validate = args

    # Get Uniprot_ID
    try:
        uni_id = id_from_ddg_path(file)
    except ValueError as e:
        logger.warning(f"Skipping ΔΔG file {file}: {e}")
        return "bad_filename"

    # Get paths of all fragments for this protein. Restrict to .csv so we don't
    # pick up sibling files (e.g. RaSP's prism_cavity_*.txt) that happen to share
    # the UniProt accession in their name.
    lst_path_prot = glob.glob(os.path.join(path_dir, f"*{uni_id}*.csv"))
    frag = len(lst_path_prot) > 1

    # Pre-validate the canonical sequence handle once (cheap). These are
    # demoted to DEBUG because, for the public RaSP bundle, they're an
    # expected consequence of UniProt-snapshot drift between the bundle and
    # the dataset — predictions for proteins Oncodrive3D doesn't analyse.
    if validate:
        if canonical_seq is None:
            logger.debug(
                f"Skipping ΔΔG for {uni_id}: UniProt ID not found in "
                "datasets/seq_for_mut_prob.tsv"
            )
            return "no_uniprot_match"
        if not isinstance(canonical_seq, str) or len(canonical_seq) == 0:
            logger.debug(
                f"Skipping ΔΔG for {uni_id}: sequence missing in "
                "datasets/seq_for_mut_prob.tsv for this UniProt ID"
            )
            return "no_seq_available"
        seq_len = len(canonical_seq)
        # Pre-convert canonical to a numpy array of single-char codes so we can
        # index it with the int-pos array per fragment.
        canonical_arr = np.frombuffer(canonical_seq.upper().encode("ascii"), dtype="S1")

    # Accumulate variants across fragments. Validation runs inline.
    ddg_dict = {}
    n_total = 0
    n_mismatch = 0
    for path_prot in lst_path_prot:
        try:
            wt_arr, pos_arr, alt_arr, ddg_arr = _parse_ddg_csv(path_prot)
        except Exception as e:
            logger.warning(
                f"Skipping ΔΔG for {uni_id}: unreadable CSV "
                f"({os.path.basename(path_prot)}): "
                f"{type(e).__name__}: {str(e)[:120]}"
            )
            return "unreadable"

        if len(pos_arr) == 0:
            # No well-formed variants in this fragment — skip it; other fragments
            # may still contribute.
            continue

        if validate:
            pos_int = pos_arr.astype(int)
            oor_mask = (pos_int < 1) | (pos_int > seq_len)
            if oor_mask.any():
                bad_pos = int(pos_int[oor_mask][0])
                logger.warning(
                    f"Skipping ΔΔG for {uni_id}: "
                    f"position {bad_pos} out of range (canonical length {seq_len})"
                )
                return "out_of_range"
            # Vectorised WT-mismatch tally. ``np.char.upper`` operates on ``S1``
            # dtype directly, so we skip the unicode round-trip.
            wt_upper = np.char.upper(wt_arr.astype("S1"))
            canonical_at_pos = canonical_arr[pos_int - 1]
            n_total += pos_int.size
            n_mismatch += int(np.count_nonzero(canonical_at_pos != wt_upper))

        # Accumulate into the position → {alt → ddg} dict. The inner loop is
        # pure Python (numpy → dict) but no longer involves per-row regex.
        if frag:
            for p, a, d in zip(pos_arr, alt_arr, ddg_arr):
                bucket = ddg_dict.setdefault(p, {})
                bucket.setdefault(a, []).append(float(d))
        else:
            for p, a, d in zip(pos_arr, alt_arr, ddg_arr):
                ddg_dict.setdefault(p, {})[a] = float(d)

    # Final validation gates
    if validate:
        if n_total == 0:
            logger.warning(f"Skipping ΔΔG for {uni_id}: no valid variants in CSV")
            return "no_variants"
        mismatch_rate = n_mismatch / n_total
        if mismatch_rate > wt_mismatch_threshold:
            logger.warning(
                f"Skipping ΔΔG for {uni_id}: "
                f"WT mismatch rate {mismatch_rate:.1%} exceeds threshold "
                f"{wt_mismatch_threshold:.1%}"
            )
            return "wt_mismatch"

    # Average across fragments for the variants seen in more than one
    if frag:
        for pos in ddg_dict:
            for alt in ddg_dict[pos]:
                ddg_dict[pos][alt] = float(np.mean(ddg_dict[pos][alt]))

    # Save dict
    save_json(output_path, uni_id, ddg_dict)
    return None


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
    # Collect (file, uni_id), deduplicating by uni_id so fragmented proteins
    # don't spawn duplicate workers that race to write the same JSON. The
    # worker globs all fragments for its uni_id internally, so a single
    # entry per accession is sufficient.
    files_processed = set(glob.glob(os.path.join(output_path, "*.json")))
    lst_files = []
    seen_uni_ids = set()
    for file in os.listdir(input_path):
        if not file.endswith(".csv"):
            continue
        try:
            uni_id = id_from_ddg_path(file)
        except ValueError as e:
            logger.warning(f"Skipping {file}: {e}")
            continue
        if uni_id in seen_uni_ids:
            continue
        if os.path.join(output_path, f"{uni_id}_ddg.json") in files_processed:
            continue
        seen_uni_ids.add(uni_id)
        lst_files.append((file, uni_id))
    ## Save dict for each proteins
    logger.debug(f"Input: {input_path}")
    logger.debug(f"Output: {output_path}")
    if len(lst_files) > 0:
        n_proteins = len(lst_files)
        logger.info(f"Parsing ΔΔG of {n_proteins} proteins...")

        # Pass only this protein's canonical sequence to each task (avoids
        # pickling the full seq_map per worker invocation).
        args_list = [
            (file, input_path, output_path,
             seq_map.get(uni_id), wt_mismatch_threshold, validate)
            for file, uni_id in lst_files
        ]
        # imap_unordered streams results back as workers finish; we use it
        # for parent-side progress reporting. chunksize amortises pickle/
        # dispatch overhead — small enough to keep the progress signal
        # responsive on small inputs, large enough to matter on big ones.
        chunksize = max(1, min(64, n_proteins // (4 * max(1, threads))))
        report_every = max(1, n_proteins // 20)  # ~20 progress lines total
        n_done = 0
        skip_counts = Counter()
        with Pool(processes=threads) as pool:
            for result in pool.imap_unordered(parse_ddg_rasp_worker, args_list,
                                              chunksize=chunksize):
                n_done += 1
                if result is not None:
                    skip_counts[result] += 1
                if n_done % report_every == 0 or n_done == n_proteins:
                    logger.info(f"Parsed {n_done}/{n_proteins} proteins")

        # Coverage and skip-reason summary
        n_skipped = sum(skip_counts.values())
        n_produced = n_done - n_skipped
        if validate and len(seq_map) > 0:
            pct = 100.0 * n_produced / len(seq_map)
            logger.info(
                f"ΔΔG coverage: {n_produced:,} / {len(seq_map):,} proteins "
                f"in datasets/seq_for_mut_prob.tsv ({pct:.1f}%)"
            )
        else:
            logger.info(f"ΔΔG produced: {n_produced:,} JSONs")
        if skip_counts:
            logger.info(f"Skipped {n_skipped:,} proteins:")
            for key, count in skip_counts.most_common():
                label = _SKIP_REASONS.get(key, key)
                logger.info(f"  - {label}: {count:,}")
    else:
        logger.debug("DDG not found: Skipping...")

    # Remove the original folder
    logger.debug(f"Deleting {input_path}")
    rm_dir(input_path)
    logger.info("Parsing of DDG completed!")