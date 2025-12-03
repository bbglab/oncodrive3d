#!/usr/bin/env python3
"""
MANE maintenance CLI for Oncodrive3D.

Example:
    python tools/update_samplesheet_and_structures.py \
        --samplesheet-folder data/251202 \
        --predicted-dir /path/to/nfcore/pdbs \
        --canonical-dir /path/to/af_canonical_pdbs
"""
from __future__ import annotations

import gzip
import os
import socket
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import click
import pandas as pd
import yaml
from Bio import SeqIO

try:
    _AVAILABLE_CORES = len(os.sched_getaffinity(0))
except (AttributeError, OSError):
    _AVAILABLE_CORES = os.cpu_count() or 1

AVAILABLE_WORKERS = max(1, _AVAILABLE_CORES)

@dataclass
class Settings:
    """Runtime knobs controlling canonical reuse, parallelism, and filtering."""
    enable_canonical_reuse: bool
    max_workers: int
    filter_long_sequences: bool
    max_sequence_length: int


@dataclass
class PipelinePaths:
    """Resolved filesystem locations for all inputs/outputs of the workflow."""
    samplesheet_path: Path
    fasta_dir: Path
    predicted_bundle_dir: Path
    predicted_pdb_dir: Path
    missing_dir: Path
    missing_fasta_dir: Path
    missing_samplesheet_path: Path
    retrieved_dir: Path
    retrieved_pdb_dir: Path
    final_bundle_dir: Path
    final_pdb_dir: Path
    canonical_pdb_dir: Optional[Path]
    mane_missing_path: Optional[Path]
    mane_summary_path: Path
    cgc_list_path: Optional[Path]


def load_config(config_path: Path) -> dict:
    """Load the project configuration from YAML."""
    with config_path.open() as fh:
        return yaml.safe_load(fh)


def detect_environment(config: dict) -> tuple[str, dict]:
    """Detect whether we are running on the cluster or locally and select paths."""
    hostname = socket.gethostname()
    irb_prefix = config["workspace"].get("irb_hostname_prefix", "irb")
    bbg_prefix = config["workspace"].get("bbg_hostname_prefix", "bbg")
    if hostname.startswith(irb_prefix):
        env_key = "cluster"
        workspace_key = "workspace_root_irb"
    elif hostname.startswith(bbg_prefix):
        env_key = "cluster"
        workspace_key = "workspace_root_bbg"
    else:
        env_key = "local"
        workspace_key = None
    path_config = config["paths"][env_key].copy()
    path_config["_workspace_key"] = workspace_key
    return env_key, path_config


def build_paths(
    config: dict,
    config_path: Path,
    path_config: dict,
    samplesheet_folder: str,
    canonical_dir: Optional[Path],
) -> PipelinePaths:
    """Resolve all filesystem paths for the selected environment and samplesheet."""
    repo_root = config_path.parent
    workspace_key = path_config.pop("_workspace_key")
    workspace_root: Path
    if workspace_key:
        workspace_value = path_config.get(workspace_key)
        if not workspace_value:
            raise KeyError(f"Missing workspace path for key: {workspace_key}")
        workspace_root = Path(workspace_value)
    else:
        workspace_root = Path(path_config.get("workspace_root", repo_root))

    format_kwargs = {
        "samplesheet_folder": samplesheet_folder,
        "WORKSPACE": str(workspace_root),
    }

    def resolve(value: str | Path | None) -> Optional[Path]:
        """Format a template path and make it absolute relative to workspace_root."""
        if value is None:
            return None
        if isinstance(value, Path):
            path_str = str(value)
        else:
            path_str = value
        resolved = path_str.format(**format_kwargs)
        path = Path(resolved)
        if not path.is_absolute():
            path = workspace_root / path
        return path

    samplesheet_path = resolve(path_config["samplesheet_relpath"])
    fasta_dir = resolve(path_config["fasta_relpath"])
    predicted_bundle_dir = resolve(path_config["predicted_relpath"])
    missing_dir = resolve(path_config["missing_relpath"])
    retrieved_dir = resolve(path_config["retrieved_relpath"])
    final_bundle_dir = resolve(path_config["final_bundle_relpath"])
    canonical_pdb_dir = None
    if canonical_dir:
        canonical_pdb_dir = canonical_dir.resolve()
        if not canonical_pdb_dir.exists():
            raise FileNotFoundError(f"--canonical-dir not found: {canonical_pdb_dir}")
    mane_missing_path = resolve(path_config.get("mane_missing_path"))
    mane_summary_path = resolve(path_config.get("mane_summary_path"))
    cgc_list_path = resolve(path_config.get("cgc_list_path"))

    if not canonical_dir:
        canonical_pdb_dir = None

    required = {
        "samplesheet": samplesheet_path,
        "fasta_dir": fasta_dir,
        "mane_summary_path": mane_summary_path,
    }
    for label, path in required.items():
        if path is None or not path.exists():
            raise FileNotFoundError(f"Required path for {label} not found: {path}")
    if cgc_list_path and not cgc_list_path.exists():
        raise FileNotFoundError(f"Required path for cgc_list_path not found: {cgc_list_path}")

    predicted_pdb_dir = predicted_bundle_dir / "pdbs"
    missing_fasta_dir = missing_dir / "fasta"
    missing_samplesheet_path = missing_dir / "samplesheet.csv"
    retrieved_pdb_dir = retrieved_dir / "pdbs"
    final_pdb_dir = final_bundle_dir / "pdbs"

    for path in [
        predicted_bundle_dir,
        predicted_pdb_dir,
        missing_dir,
        missing_fasta_dir,
        retrieved_pdb_dir,
        final_pdb_dir,
    ]:
        path.mkdir(parents=True, exist_ok=True)

    paths = PipelinePaths(
        samplesheet_path=samplesheet_path,
        fasta_dir=fasta_dir,
        predicted_bundle_dir=predicted_bundle_dir,
        predicted_pdb_dir=predicted_pdb_dir,
        missing_dir=missing_dir,
        missing_fasta_dir=missing_fasta_dir,
        missing_samplesheet_path=missing_samplesheet_path,
        retrieved_dir=retrieved_dir,
        retrieved_pdb_dir=retrieved_pdb_dir,
        final_bundle_dir=final_bundle_dir,
        final_pdb_dir=final_pdb_dir,
        canonical_pdb_dir=canonical_pdb_dir,
        mane_missing_path=mane_missing_path,
        mane_summary_path=mane_summary_path,
        cgc_list_path=cgc_list_path,
    )
    return paths


def ensure_dir(path: Path) -> Path:
    """Create a directory (and parents) if it does not exist."""
    path.mkdir(parents=True, exist_ok=True)
    return path


def list_predicted_sequences(pdb_dir: Path) -> set[str]:
    """Return the ENSP IDs inferred from the filenames inside pdb_dir."""
    if not pdb_dir.exists():
        return set()
    return {p.name.split(".")[0] for p in pdb_dir.glob("*.pdb*")}


def copy_fastas(df: pd.DataFrame, src_dir: Path, dst_dir: Path) -> pd.DataFrame:
    """Copy FASTA files for the requested ENSPs into dst_dir and return their paths."""
    ensure_dir(dst_dir)
    records = []
    for ens_id in df["sequence"]:
        src = src_dir / f"{ens_id}.fasta"
        if not src.exists():
            print(f"[WARNING] Missing FASTA for {ens_id}: {src}")
            continue
        dst = dst_dir / src.name
        shutil.copy2(src, dst)
        records.append({"sequence": ens_id, "fasta": dst.as_posix()})
    return pd.DataFrame(records)


def extract_seq_from_pdb(pdb_path: Path) -> str:
    """Extract the first sequence from a PDB/PDB.GZ file."""
    opener = gzip.open if pdb_path.suffix == ".gz" else open
    with opener(pdb_path, "rt") as handle:
        seq_records = [str(record.seq) for record in SeqIO.parse(handle, "pdb-seqres")]
    return seq_records[0] if seq_records else ""


def index_canonical_pdbs(pdb_dir: Path, max_workers: int) -> pd.DataFrame:
    """Build a dataframe describing every canonical MANE structure in pdb_dir."""
    files = [entry for entry in pdb_dir.iterdir() if entry.is_file() and entry.name.startswith("AF-")]
    rows = []
    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        futures = {pool.submit(extract_seq_from_pdb, path): path for path in files}
        for future in as_completed(futures):
            path = futures[future]
            seq = future.result()
            parts = path.name.split("-")
            if len(parts) < 3:
                continue
            uniprot_id = parts[1]
            fragment_token = parts[2]
            fragment = fragment_token.replace("F", "")
            rows.append(
                {
                    "path": path.as_posix(),
                    "uniprot_id": uniprot_id,
                    "fragment": fragment,
                    "seq_len": len(seq),
                }
            )
    return pd.DataFrame(rows)


def merge_structure_bundles(bundle_dirs: Iterable[Path], output_dir: Path) -> pd.DataFrame:
    """Merge multiple ENSP bundles into a single `output_dir` with deduplicated samplesheet."""
    output_dir = ensure_dir(output_dir)
    target_pdb_dir = ensure_dir(output_dir / "pdbs")
    merged = []
    for bundle_dir in bundle_dirs:
        samplesheet_file = bundle_dir / "samplesheet.csv"
        pdb_dir = bundle_dir / "pdbs"
        if not samplesheet_file.exists() or not pdb_dir.exists():
            print(f"[SKIP] {bundle_dir} is missing samplesheet or pdbs/")
            continue
        df = pd.read_csv(samplesheet_file)
        merged.append(df)
        for pdb_file in pdb_dir.glob("*.pdb*"):
            shutil.copy2(pdb_file, target_pdb_dir / pdb_file.name)
    if not merged:
        raise RuntimeError("No valid bundles provided for merging.")
    combined = pd.concat(merged, ignore_index=True).drop_duplicates(subset=["sequence"])
    combined.to_csv(output_dir / "samplesheet.csv", index=False)
    print(f"Merged {len(merged)} bundles → {len(combined)} unique ENSP entries")
    return combined


def sync_predicted_bundle(src_dir: Path, bundle_dir: Path) -> pd.DataFrame:
    """Copy nf-core outputs from src_dir into the curated bundle layout and emit samplesheet."""
    src_dir = Path(src_dir)
    bundle_dir = ensure_dir(bundle_dir)
    if not src_dir.exists():
        raise FileNotFoundError(f"predicted_dir not found: {src_dir}")
    if src_dir.resolve() == bundle_dir.resolve():
        print(f"[INFO] predicted_dir already points to {bundle_dir}. Skipping sync.")
        samplesheet_path = bundle_dir / "samplesheet.csv"
        if samplesheet_path.exists():
            return pd.read_csv(samplesheet_path)
        return pd.DataFrame()

    dst_pdb_dir = ensure_dir(bundle_dir / "pdbs")
    records = []
    for pdb_file in sorted(src_dir.glob("*.pdb*")):
        ens_id = pdb_file.name.split(".")[0]
        shutil.copy2(pdb_file, dst_pdb_dir / pdb_file.name)
        records.append({"sequence": ens_id})

    if not records:
        print(f"[WARNING] No predicted PDBs found in {src_dir}")
        return pd.DataFrame()

    df = pd.DataFrame(records).drop_duplicates().sort_values("sequence")
    df.to_csv(bundle_dir / "samplesheet.csv", index=False)
    print(f"Synced {len(df)} predicted structures → {bundle_dir}")
    return df


def prune_missing_bundle(missing_dir: Path, sequences: set[str]) -> None:
    """Remove already satisfied ENSPs from the missing samplesheet/fasta directories."""
    if not sequences:
        return
    samplesheet_file = missing_dir / "samplesheet.csv"
    if not samplesheet_file.exists():
        print(f"[INFO] Missing samplesheet for pruning: {samplesheet_file}")
        return
    df = pd.read_csv(samplesheet_file)
    mask = df["sequence"].isin(sequences)
    if not mask.any():
        return
    kept = df[~mask].copy()
    kept.to_csv(samplesheet_file, index=False)

    fasta_dir = missing_dir / "fasta"
    for seq in sequences:
        fasta_path = fasta_dir / f"{seq}.fasta"
        if fasta_path.exists():
            fasta_path.unlink()

    print(f"Pruned {mask.sum()} entries from missing bundle after canonical reuse")


def load_cgc_symbols(cgc_path: Optional[Path]) -> set[str]:
    """Return the set of CGC gene symbols (and synonyms) if the file exists."""
    if not cgc_path or not cgc_path.exists():
        print("[INFO] CGC list not provided; CGC prioritisation disabled.")
        return set()
    df = pd.read_csv(cgc_path, sep="\t")
    symbols = set(df.get("Gene Symbol", pd.Series(dtype=str)).dropna().str.strip())
    if "Synonyms" in df.columns:
        synonyms = df["Synonyms"].dropna().str.split(",").explode().str.strip()
        symbols.update(sym for sym in synonyms if sym)
    return {sym for sym in symbols if sym}


def compute_fasta_lengths(fasta_paths: pd.Series) -> pd.Series:
    """Compute sequence lengths for each FASTA file path."""
    lengths = []
    for path_str in fasta_paths.fillna(""):
        path = Path(path_str)
        if not path.exists():
            print(f"[WARNING] FASTA file missing for {path}")
            lengths.append(pd.NA)
            continue
        length = 0
        with path.open() as handle:
            for line in handle:
                if line.startswith(">"):
                    continue
                length += len(line.strip())
        lengths.append(length)
    return pd.Series(lengths, index=fasta_paths.index, dtype="Int64")


def filter_long_sequences(df: pd.DataFrame, missing_dir: Path, max_len: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Drop ENSPs whose FASTA length exceeds max_len and log the removed entries."""
    if "length" not in df.columns:
        return df, pd.DataFrame()
    mask = df["length"].fillna(-1) > max_len
    if not mask.any():
        return df, pd.DataFrame()
    kept = df[~mask].copy()
    removed = df[mask].copy()

    removed_clean = removed.drop(columns=["symbol", "CGC", "length"], errors="ignore")
    removed_path = missing_dir / "samplesheet_removed_long.csv"
    removed_clean.to_csv(removed_path, index=False)

    fasta_dir = missing_dir / "fasta"
    for seq in removed["sequence"].unique():
        fasta_path = fasta_dir / f"{seq}.fasta"
        if fasta_path.exists():
            fasta_path.unlink()
        else:
            print(f"[WARNING] FASTA file to remove not found: {fasta_path}")

    print(f"Removed {len(removed_clean)} entries longer than {max_len} residues")
    return kept, removed_clean


def annotate_missing_with_cgc(
    missing_sheet: Path,
    cgc_path: Optional[Path],
    mane_summary_path: Path,
) -> pd.DataFrame:
    """Annotate the missing samplesheet with CGC tags and computed lengths."""
    missing_df = pd.read_csv(missing_sheet)
    if missing_df.empty:
        print("[INFO] No entries left in missing samplesheet.")
        return missing_df

    cgc_symbols = load_cgc_symbols(cgc_path)
    mane_summary = pd.read_table(mane_summary_path)

    column_aliases = {
        "ensembl_prot": {"Ensembl_prot", "ensembl_prot", "Ens_Prot_ID"},
        "refseq_prot": {"RefSeq_prot", "refseq_prot"},
        "symbol": {"symbol", "Gene Symbol", "gene_symbol"},
    }

    rename_map = {}
    for target, candidates in column_aliases.items():
        for candidate in candidates:
            if candidate in mane_summary.columns:
                rename_map[candidate] = target
                break
    required = {"ensembl_prot", "refseq_prot", "symbol"}
    if not required.issubset(rename_map.values()):
        missing = required - set(rename_map.values())
        raise KeyError(f"Missing columns in MANE summary: {missing}")

    mane_summary = mane_summary.rename(columns=rename_map)
    annotations = mane_summary[["ensembl_prot", "refseq_prot", "symbol"]].drop_duplicates()
    annotations["CGC"] = annotations["symbol"].isin(cgc_symbols).astype(int)

    seq_map = annotations.dropna(subset=["ensembl_prot"]).set_index("ensembl_prot")
    refseq_map = annotations.dropna(subset=["refseq_prot"]).set_index("refseq_prot")

    missing_df["symbol"] = missing_df["sequence"].map(seq_map["symbol"] if not seq_map.empty else pd.Series(dtype=str))
    if "refseq_prot" in missing_df.columns:
        missing_df["symbol"] = missing_df["symbol"].fillna(
            missing_df["refseq_prot"].map(refseq_map["symbol"] if not refseq_map.empty else pd.Series(dtype=str))
        )

    missing_df["CGC"] = missing_df["sequence"].map(seq_map["CGC"] if not seq_map.empty else pd.Series(dtype=int))
    if "refseq_prot" in missing_df.columns:
        missing_df["CGC"] = missing_df["CGC"].fillna(
            missing_df["refseq_prot"].map(refseq_map["CGC"] if not refseq_map.empty else pd.Series(dtype=int))
        )

    missing_df["CGC"] = missing_df["CGC"].fillna(0).astype(int)
    missing_df["symbol"] = missing_df["symbol"].fillna("")

    if "fasta" in missing_df.columns:
        missing_df["length"] = compute_fasta_lengths(missing_df["fasta"])
    else:
        missing_df["length"] = pd.NA

    sort_columns = ["length", "sequence"]
    ascending = [True, True]
    if cgc_symbols:
        sort_columns = ["CGC"] + sort_columns
        ascending = [False] + ascending

    missing_df = missing_df.sort_values(sort_columns, ascending=ascending).reset_index(drop=True)

    return missing_df


def reuse_canonical_structures(
    samplesheet_missing: pd.DataFrame,
    paths: PipelinePaths,
    settings: Settings,
) -> pd.DataFrame:
    """Copy canonical MANE structures for matching ENSPs and prune them from the missing set."""
    if not settings.enable_canonical_reuse:
        print("Canonical reuse disabled by config.")
        return pd.DataFrame()
    if not paths.canonical_pdb_dir or not paths.canonical_pdb_dir.exists():
        print("Canonical PDB directory not found. Skipping reuse.")
        return pd.DataFrame()
    if not paths.mane_missing_path or not paths.mane_missing_path.exists():
        print("mane_missing.csv not found. Skipping reuse.")
        return pd.DataFrame()

    if "refseq_prot" not in samplesheet_missing.columns:
        print("[INFO] Missing refseq_prot column; cannot perform canonical reuse.")
        return pd.DataFrame()

    canonical_index = index_canonical_pdbs(paths.canonical_pdb_dir, settings.max_workers)
    if canonical_index.empty:
        print(
            f"[ERROR] No canonical PDB files found in {paths.canonical_pdb_dir}. "
            "Verify --canonical-dir points to the AlphaFold DB download."
        )
        return pd.DataFrame()
    print(f"Indexed {len(canonical_index):,} canonical PDB files")

    mane_missing_df = pd.read_csv(paths.mane_missing_path)
    rename_map = {}
    if "uniprot_accession(s)" in mane_missing_df.columns:
        rename_map["uniprot_accession(s)"] = "uniprot_id"
    if "RefSeq_prot" in mane_missing_df.columns:
        rename_map["RefSeq_prot"] = "refseq_prot"
    mane_missing_df = mane_missing_df.rename(columns=rename_map)
    if not {"uniprot_id", "refseq_prot"}.issubset(mane_missing_df.columns):
        print("[INFO] Cannot map UniProt↔RefSeq from mane_missing.csv. Skipping reuse.")
        return pd.DataFrame()
    mane_missing_df = mane_missing_df[["uniprot_id", "refseq_prot"]].dropna().drop_duplicates()

    canonical_with_refseq = canonical_index.merge(mane_missing_df, on="uniprot_id", how="left").dropna(subset=["refseq_prot"])
    retrievable = samplesheet_missing.merge(canonical_with_refseq, on="refseq_prot", how="inner")
    print(f"Canonical structures available for {len(retrievable):,} ENSP entries")

    retrieved_records = []
    ensure_dir(paths.retrieved_pdb_dir)
    for _, row in retrievable.iterrows():
        src_file = Path(row["path"])
        fragment = row.get("fragment") or "1"
        dst_name = f"{row.sequence}.{fragment}.alphafold.pdb"
        dst_file = paths.retrieved_pdb_dir / dst_name
        try:
            if src_file.suffix == ".gz":
                with gzip.open(src_file, "rt") as src, open(dst_file, "wt") as dst:
                    shutil.copyfileobj(src, dst)
            else:
                shutil.copy2(src_file, dst_file)
            retrieved_records.append({"sequence": row.sequence})
        except Exception as exc:  # pylint: disable=broad-except
            print(f"[WARNING] Failed to process {src_file} → {dst_file}: {exc}")

    if not retrieved_records:
        print("[INFO] No canonical structures retrieved.")
        return pd.DataFrame()

    retrieved_df = pd.DataFrame(retrieved_records).drop_duplicates()
    retrieved_df.to_csv(paths.retrieved_dir / "samplesheet.csv", index=False)
    print(f"Stored {len(retrieved_df):,} canonical structures in {paths.retrieved_pdb_dir}")
    prune_missing_bundle(paths.missing_dir, set(retrieved_df["sequence"]))
    return retrieved_df


def run_pipeline(
    paths: PipelinePaths,
    settings: Settings,
    predicted_raw_dir: Optional[Path],
    canonical_dir: Optional[Path],
) -> None:
    """Execute the MANE maintenance workflow end-to-end."""
    samplesheet = pd.read_csv(paths.samplesheet_path)
    print(f"Loaded {len(samplesheet):,} samples from {paths.samplesheet_path}")

    if predicted_raw_dir:
        if not predicted_raw_dir.exists():
            raise FileNotFoundError(f"--predicted-dir not found: {predicted_raw_dir}")
        sync_predicted_bundle(predicted_raw_dir, paths.predicted_bundle_dir)
    else:
        print("[INFO] No --predicted-dir supplied; skipping nf-core sync and using existing predicted bundle.")

    predicted_ids = list_predicted_sequences(paths.predicted_pdb_dir)
    print(f"Found {len(predicted_ids):,} predicted ENSP structures in {paths.predicted_pdb_dir}")

    samplesheet_pred = samplesheet[samplesheet.sequence.isin(predicted_ids)].copy()
    samplesheet_missing = samplesheet[~samplesheet.sequence.isin(predicted_ids)].copy()
    print(f"Already predicted: {len(samplesheet_pred):,}")
    print(f"Still missing:    {len(samplesheet_missing):,}")

    missing_fastas_df = copy_fastas(samplesheet_missing, paths.fasta_dir, paths.missing_fasta_dir)
    missing_fastas_df.to_csv(paths.missing_samplesheet_path, index=False)
    print(f"Prepared {len(missing_fastas_df):,} FASTA files for nf-core at {paths.missing_samplesheet_path}")
    if not missing_fastas_df.empty:
        samplesheet_missing = samplesheet_missing[samplesheet_missing.sequence.isin(missing_fastas_df["sequence"])]
    else:
        samplesheet_missing = samplesheet_missing.iloc[0:0].copy()

    retrieved_df = pd.DataFrame()
    if settings.enable_canonical_reuse:
        retrieved_df = reuse_canonical_structures(samplesheet_missing, paths, settings)
        if retrieved_df.empty:
            print("[INFO] Skipping canonical reuse because no PDBs were harvested.")
    else:
        print("[INFO] No canonical directory provided; skipping canonical reuse.")

    annotated_missing_df = annotate_missing_with_cgc(
        paths.missing_samplesheet_path,
        paths.cgc_list_path,
        paths.mane_summary_path,
    )

    removed_long_df = pd.DataFrame()
    if settings.filter_long_sequences and not annotated_missing_df.empty:
        annotated_missing_df, removed_long_df = filter_long_sequences(
            annotated_missing_df, paths.missing_dir, settings.max_sequence_length
        )

    clean_missing_df = annotated_missing_df.drop(columns=["symbol", "CGC", "length"], errors="ignore")
    clean_missing_df.to_csv(paths.missing_samplesheet_path, index=False)
    print(f"Updated missing samplesheet saved to {paths.missing_samplesheet_path}")

    if settings.filter_long_sequences and not removed_long_df.empty:
        print(f"Long sequences recorded at {paths.missing_dir / 'samplesheet_removed_long.csv'}")

    bundles_to_merge = []
    if paths.predicted_pdb_dir.exists() and any(paths.predicted_pdb_dir.glob("*.pdb*")):
        bundles_to_merge.append(paths.predicted_bundle_dir)
    if paths.retrieved_pdb_dir.exists() and any(paths.retrieved_pdb_dir.glob("*.pdb*")):
        bundles_to_merge.append(paths.retrieved_dir)

    if not bundles_to_merge:
        print("[INFO] No predicted/retrieved bundles available; skipping final merge.")
        return

    merge_structure_bundles(bundles_to_merge, paths.final_bundle_dir)
    print(f"Final bundle written to {paths.final_bundle_dir}")


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--samplesheet-folder", required=True, help="Folder containing the MANE samplesheet/fasta bundle created by Oncodrive3D.")
@click.option("--config-path", default="config.yaml", type=click.Path(path_type=Path), show_default=True)
@click.option(
    "--predicted-dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Directory with nf-core/proteinfold predicted ENSP PDBs to sync (optional).",
)
@click.option(
    "--canonical-dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Directory containing the AlphaFold DB canonical PDBs (optional, enables reuse).",
)
@click.option(
    "--max-workers",
    type=click.IntRange(1, AVAILABLE_WORKERS),
    default=AVAILABLE_WORKERS,
    show_default=True,
    help="Parallel workers.",
)
@click.option(
    "--filter-long-sequences/--no-filter-long-sequences",
    default=True,
    show_default=True,
    help="Remove entries whose FASTA length exceeds --max-sequence-length from the nf-core input.",
)
@click.option(
    "--max-sequence-length",
    type=int,
    default=2700,
    show_default=True,
    help="Maximum sequence length (residues) allowed in the nf-core missing set.",
)
def cli(
    samplesheet_folder: str,
    config_path: Path,
    predicted_dir: Optional[Path],
    canonical_dir: Optional[Path],
    max_workers: int,
    filter_long_sequences: bool,
    max_sequence_length: int,
) -> None:
    """CLI wrapper around the MANE maintenance workflow."""
    if not predicted_dir and not canonical_dir:
        raise click.UsageError("Provide at least --predicted-dir or --canonical-dir so there is work to perform.")

    config = load_config(config_path)
    _, path_config = detect_environment(config)
    paths = build_paths(
        config,
        config_path,
        path_config,
        samplesheet_folder,
        canonical_dir=canonical_dir,
    )
    settings = Settings(
        enable_canonical_reuse=bool(paths.canonical_pdb_dir),
        max_workers=max_workers,
        filter_long_sequences=filter_long_sequences,
        max_sequence_length=max_sequence_length,
    )
    run_pipeline(paths, settings, predicted_dir, canonical_dir)


if __name__ == "__main__":
    cli()
