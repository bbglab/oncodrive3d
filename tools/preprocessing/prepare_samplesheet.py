"""
prepare_samplesheet.py

Utility to download the MANE Ensembl protein FASTA, filter out entries
already in AlphaFold, write individual .fasta files, and assemble a
samplesheet.csv for downstream nf-core AlphaFold pipeline.

Example usage:
    python -m tools.preprocessing.prepare_samplesheet \
        --datasets-dir  /data/bbg/nobackup/scratch/oncodrive3d/datasets_mane_240506/ \
        --output-dir    /data/bbg/nobackup/scratch/oncodrive3d/mane_missing/data
"""


import os
import click
import gzip
import time
import pandas as pd
from pathlib import Path
from scripts.datasets.utils import download_single_file
# import logging

# logger = logging.getLogger(__name__)

class ManeSamplesheetBuilder:
    """
    Download MANE FASTA, merge with AlphaFold mapping, 
    write out individual FASTAs, and assemble the final samplesheet.
    """

    def __init__(
        self,
        datasets_dir: str,
        output_dir: str,
        mane_version: str = "1.4",
        max_attempts: int = 15,
        no_fragments: bool = True,
        cores: int = 1,
    ):
        self.datasets_dir = Path(datasets_dir)
        self.output_dir = Path(output_dir)
        self.mane_version = mane_version
        self.max_attempts = max_attempts
        self.no_fragments = no_fragments
        self.cores = cores
        
        self.fasta_gz = self.output_dir / f"MANE.GRCh38.v{mane_version}.ensembl_protein.faa.gz"
        self.fasta_dir = self.output_dir / "fasta"
        self.samplesheet_path = self.output_dir / "samplesheet.csv"

        self.output_dir.mkdir(parents=True, exist_ok=True)


    def _download_mane_fasta(self):
        """Download the MANE Ensembl protein FASTA (gzipped)."""
        url = (
            f"https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human"
            f"/release_{self.mane_version}/MANE.GRCh38.v{self.mane_version}.ensembl_protein.faa.gz"
        )
        attempts = 0
        while attempts < self.max_attempts:
            try:
                download_single_file(url, str(self.fasta_gz), threads=self.cores)
                print("Downloaded MANE FASTA successfully.")
                return
            except Exception as e:
                attempts += 1
                print(f"Attempt {attempts} failed: {e!r}")
                time.sleep(5)
        raise RuntimeError(f"Failed to download MANE FASTA after {self.max_attempts} attempts.")


    def _parse_ncbi_mane_fasta(self) -> pd.DataFrame:
        """
        Parse the gzipped FASTA into a DataFrame with columns
        ['sequence', 'refseq'] where 'sequence' is the Ensembl ID
        (no version) and 'refseq' is the AA string.
        """
        if not self.fasta_gz.exists():
            raise FileNotFoundError(self.fasta_gz)
        ids, seqs = [], []
        with gzip.open(self.fasta_gz, "rt") as fh:
            header = None
            seq_parts = []
            for line in fh:
                line = line.rstrip()
                if line.startswith(">"):
                    if header is not None:
                        ids.append(header)
                        seqs.append("".join(seq_parts))
                    raw = line[1:].split()[0]
                    header = raw.split(".", 1)[0]
                    seq_parts = []
                else:
                    seq_parts.append(line)
            # last record
            if header is not None:
                ids.append(header)
                seqs.append("".join(seq_parts))

        return pd.DataFrame({"sequence": ids, "refseq": seqs})


    def _load_mane_af(self) -> pd.DataFrame:
        """
        Load the mapping of RefSeq-UniProt-Ensembl_prot
        and return a DataFrame with 'refseq_prot' and
        version‐stripped 'Ensembl_prot'.
        """
        map_csv = self.datasets_dir / "mane_refseq_prot_to_alphafold.csv"
        sum_gz   = self.datasets_dir / "mane_summary.txt.gz"

        df_map = pd.read_csv(map_csv)
        self.df_mane_summary = (
            pd.read_csv(sum_gz, compression="gzip", sep="\t")
            .rename(columns={"RefSeq_prot": "refseq_prot"})
            .dropna(subset=["Ensembl_prot", "refseq_prot"])
            )
        df = df_map.merge(self.df_mane_summary[["refseq_prot", "Ensembl_prot"]], on="refseq_prot", how="inner")
        df["Ensembl_prot"] = df["Ensembl_prot"].str.split(".", n=1).str[0]
        
        return df


    def _add_refseq(self, df):
        """Add RefSeq Protein information to the samplesheet"""
        mane_summary = self.df_mane_summary.rename(
            columns={"Ensembl_prot": "sequence"}
            )[["sequence", "refseq_prot"]]
        mane_summary["sequence"] = mane_summary["sequence"].str.split(".", n=1).str[0]
        df = df.merge(mane_summary, on="sequence", how="left")
        
        return df


    def write_fastas_and_update_sheet(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        For every row in df, write a .fasta file under self.fasta_dir,
        then insert 'fasta' and 'length' columns and return samplesheet df.
        """

        df = df.copy()
        df["length"] = df["refseq"].str.len()
        self.fasta_dir.mkdir(exist_ok=True, parents=True)

        if self.no_fragments:
            df = df[df["length"] <= 2700].reset_index(drop=True)

        # Build fasta paths and write files
        fasta_paths = []
        for seq_id, seq, length in zip(df["sequence"], df["refseq"], df["length"]):
            p = self.fasta_dir / f"{seq_id}.fasta"
            fasta_str = f">{seq_id} | {length} aa\n{seq}\n"
            p.write_text(fasta_str)
            fasta_paths.append(str(p))

        df.insert(1, "fasta", fasta_paths)

        return df


    def build(self) -> pd.DataFrame:
        """
        Run the full workflow and write out fastas and samplesheet.csv.
        """
        self._download_mane_fasta()
        mane_all = self._parse_ncbi_mane_fasta()
        mane_af = self._load_mane_af()

        # Keep only those not already in AF
        mane_missing_af = mane_all[~mane_all["sequence"].isin(mane_af["Ensembl_prot"])]
        sheet = self.write_fastas_and_update_sheet(mane_missing_af)
        
        # Save
        sheet = self._add_refseq(sheet)
        sheet.to_csv(self.samplesheet_path, index=False)
        print(f"Wrote samplesheet to {self.samplesheet_path}")
        
        return sheet


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--datasets-dir", "-d",
    required=True,         # This could be optional if we download mane_refseq_prot_to_alphafold.csv from AlphaFold DB and mane_summary.txt.gz from ncbi
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Path to Oncodrive3D datasets or folder containing AF‑MANE mapping (mane_refseq_prot_to_alphafold.csv) and summary (mane_summary.txt.gz) files"
)
@click.option(
    "--output-dir", "-o",
    required=True,
    type=click.Path(file_okay=False, dir_okay=True),
    help="Where to write the downloaded FASTA, individual .fasta files, and samplesheet"
)
@click.option(
    "--mane-version", "-v",
    default="1.4",
    show_default=True,
    help="MANE release version to download (e.g. 1.4)"
)
@click.option(
    "--no-fragments",
    is_flag=True,
    default=False,
    help="Drop any sequences ≥2700 aa from the samplesheet"
)
@click.option(
    "-c", "--cores",
    type=click.IntRange(min=1, max=len(os.sched_getaffinity(0))),
    default=len(os.sched_getaffinity(0)),
    show_default=True,
    help="Number of cores to use in the computation"
)
def main(datasets_dir, output_dir, mane_version, no_fragments, cores):
    """
    Build a nf‑core AlphaFold samplesheet by downloading MANE FASTA,
    filtering out existing AF entries, writing per‑protein FASTAs,
    and emitting a samplesheet.csv.
    """
    
    # Log the parameters
    print("Running with parameters:")
    for name, val in {
        "datasets_dir   ": datasets_dir,
        "output_dir     ": output_dir,
        "mane_version   ": mane_version,
        "no_fragments   ": no_fragments,
        "cores          ": cores
        }.items():
        print(f"{name} = {val}")

    builder = ManeSamplesheetBuilder(
        datasets_dir=datasets_dir,
        output_dir=output_dir,
        mane_version=mane_version,
        no_fragments=no_fragments
    )
    builder.build()
    

if __name__ == "__main__":
    main()
    
    # # For debugging
    # builder = ManeSamplesheetBuilder(
    #     datasets_dir="/data/bbg/nobackup/scratch/oncodrive3d/datasets_mane_240506/",
    #     output_dir="/data/bbg/nobackup/scratch/oncodrive3d/mane_missing/data/250728-all_proteins",
    #     cores=len(os.sched_getaffinity(0)),
    #     no_fragments=True
    #     )
    # builder.build()