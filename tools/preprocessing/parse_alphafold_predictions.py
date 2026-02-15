#!/usr/bin/env python3
"""
Parse nested AlphaFold predictions into a simple PDB folder.

A recent AlphaFold pipeline produces one directory per protein (ENSP...) 
with the actual structure nested inside and named `ranked_0.pdb`.

This script scans a `--predicted-dir` containing multiple `ENSP*` folders,
finds each folder's `ranked_0.pdb` (or `ranked_0.pdb.gz`), and copies it to
`--output-dir` as:

    <ENSP_ID>.1.alphafold.pdb
"""

from __future__ import annotations

import gzip
import shutil
from pathlib import Path

import click


def find_ranked_structure(ensp_folder: Path) -> Path | None:
    """Return the path to ranked_0.pdb (or ranked_0.pdb.gz) under ensp_folder."""
    candidates = list(ensp_folder.rglob("ranked_0.pdb"))
    if not candidates:
        candidates = list(ensp_folder.rglob("ranked_0.pdb.gz"))
    if not candidates:
        return None

    ensp_id = ensp_folder.name
    preferred = [path for path in candidates if path.parent.name == ensp_id]
    if preferred:
        candidates = preferred

    candidates.sort(key=lambda path: path.as_posix())
    return candidates[0]


def copy_pdb(src: Path, dst: Path) -> None:
    """Copy src -> dst, transparently handling gzipped input."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.suffix == ".gz":
        with gzip.open(src, "rb") as handle_src, dst.open("wb") as handle_dst:
            shutil.copyfileobj(handle_src, handle_dst)
    else:
        shutil.copy2(src, dst)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--predicted-dir",
    type=click.Path(path_type=Path, exists=True, file_okay=False),
    required=True,
    help="Directory containing nested AlphaFold predictions (top-level ENSP* folders).",
)
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path, file_okay=False),
    required=True,
    help="Directory where flattened PDBs will be written.",
)
def cli(predicted_dir: Path, output_dir: Path) -> None:
    """Collect ranked_0 PDBs from nested outputs and rename them by ENSP."""
    predicted_dir = predicted_dir.expanduser().resolve()
    output_dir = output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    ensp_folders = sorted(
        folder for folder in predicted_dir.iterdir() if folder.is_dir() and folder.name.startswith("ENSP")
    )
    if not ensp_folders:
        raise click.ClickException(f"No ENSP* folders found in: {predicted_dir}")

    copied = 0
    missing = 0
    for folder in ensp_folders:
        ensp_id = folder.name
        ranked = find_ranked_structure(folder)
        if ranked is None:
            click.echo(f"[WARNING] {ensp_id}: ranked_0.pdb not found")
            missing += 1
            continue

        dst = output_dir / f"{ensp_id}.1.alphafold.pdb"
        try:
            copy_pdb(ranked, dst)
            copied += 1
        except OSError as exc:
            click.echo(f"[WARNING] {ensp_id}: failed to copy {ranked} -> {dst}: {exc}")
            missing += 1

    click.echo(f"Scanned {len(ensp_folders)} ENSP folders.")
    click.echo(f"Copied {copied} structures into {output_dir}.")
    if missing:
        click.echo(f"Missing/failed: {missing}")


if __name__ == "__main__":
    cli()

