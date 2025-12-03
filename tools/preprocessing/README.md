# MANE Preprocessing Toolkit

These tools are provided for users that want to run Oncodrive3D using only structures associated to MANE Select transcripts while still covering as many proteins as possible. Infact, AlphaFold database MANE download bundle does not yet contain structures for every MANE Select transcript. Oncodrive3D relies on that bundle when building its datasets, so missing structures translate into genes that cannot be analyzed. The scripts in `tools/preprocessing/` close this gap:

- `prepare_samplesheet.py` scans the full MANE release and emits `samplesheet.csv` plus per-ENSP FASTAs for every MANE structure that is absent from the AlphaFold MANE download.
- `update_samplesheet_and_structures.py` removes the MANE entries already covered by the AlphaFold canonical bundle, reuses those canonical structures when available, and folds nf-core predictions into the custom bundle while pruning fulfilled entries from the next proteinfold run.

Together they allow to iteratively update the MANE structures and feed them back into `oncodrive3d build-datasets --custom_mane_pdb_dir <path/to/final_bundle/pdbs> --custom_mane_metadata_path  <path/to/final_bundle/samplesheet.csv>`.

## Prerequisites

Run `oncodrive3d build-datasets --mane_only` to generate the MANE mapping files consumed by `prepare_samplesheet.py`. If you plan to reuse canonical structures, run a separate `oncodrive3d build-datasets --mane` (or default `oncodrive3d build-datasets`) so you also have the AlphaFold canonical bundle whose `pdb_structures/` matching the missing MANE structures can be retrieved by the tools.

After running the tools and predicting the missing structures using the nf-core/proteinfold pipeline, rerun `oncodrive3d build-datasets --mane_only --custom_mane_pdb_dir ... --custom_mane_metadata_path ...` to inject the curated bundle into the final MANE-only datasets.

## Tool overview

| Script | Goal | When to run | Key outputs |
| --- | --- | --- | --- |
| `prepare_samplesheet.py` | Detect MANE ENSPs missing from the AlphaFold download and produce the FASTA + samplesheet inputs for nf-core/proteinfold. | Whenever there is the need to refresh the missing-set PDBs. | `MANE.GRCh38.vX.Y.ensembl_protein.faa.gz`, `fasta/*.fasta`, `samplesheet.csv`. |
| `update_samplesheet_and_structures.py` | Reuse canonical AlphaFold structures, sync nf-core predictions, and maintain the final custom MANE bundle + updated missing set. | After preparing the samplesheet, before and/or after nf-core runs. | `retrieved/`, `predicted/`, `final_bundle/`, refreshed `missing/` tree. |

---

## Installation
Requires:
- Python 3.10+
- Oncodrive3D datasets built as described in the prerequisites (a `--mane_only` run for the MANE-only baseline plus, if needed, a standard `--mane` or default run for the canonical bundle). These builds supply `mane_refseq_prot_to_alphafold.csv`, `mane_summary.txt.gz`, `mane_missing.csv`, and (optionally) the canonical PDB bundle.

```bash
uv sync
source .venv/bin/activate
```

## `prepare_samplesheet.py`

### Purpose
Download the requested MANE protein FASTA from NCBI and materialize FASTA files plus a `samplesheet.csv` containing only the ENSPs that are missing from the MANE AlphaFold download. The resulting folder is a drop-in input for nf-core/proteinfold. It relies on the MANE mapping files that ship with `oncodrive3d build-datasets --mane_only`.

### Usage

```bash
python -m tools.preprocessing.prepare_samplesheet \
    --datasets-dir /path/to/o3d_mane_only_dataset \
    --output-dir   /path/to/mane_missing/data
```

Key options:

- `--datasets-dir/-d` (**required**): Directory with the AlphaFold MANE mapping files produced by `oncodrive3d build-datasets`.
- `--output-dir/-o` (**required**): Destination folder holding the downloaded FASTA, the per-ENSP FASTAs, and `samplesheet.csv`.
- `--mane-version/-v`: Release to download from NCBI (default `1.4`).
- `--no-fragments`: Drop proteins longer than 2,700 aa (enabled by default when instantiating the builder; pass `--no-fragments/--no-no-fragments` to toggle).
- `--cores`: Number of download threads (defaults to all available cores).

### Workflow

1. Download `MANE.GRCh38.v<version>.ensembl_protein.faa.gz` from NCBI (with retries).
2. Parse the FASTA, strip version suffixes, and build a DataFrame of MANE proteins.
3. Load the AlphaFold MANE ↔ RefSeq mapping, drop any ENSPs already covered in the AlphaFold MANE download, and annotate the remaining sequences with RefSeq metadata.
4. Write per-protein FASTAs under `fasta/` and annotate sequence lengths and RefSeq IDs.
5. Emit the curated `samplesheet.csv` containing only the missing MANE entries, ready for nf-core/proteinfold.

### Outputs

Inside the chosen `--output-dir` you will find:

- `MANE.GRCh38.vX.Y.ensembl_protein.faa.gz` – cached MANE FASTA download.
- `fasta/ENSP....fasta` – one FASTA per missing ENSP (filtered for length if requested).
- `samplesheet.csv` – columns include `sequence`, `fasta`, `refseq`, `length`, `refseq_prot`, etc.

Feed these files to `update_samplesheet_and_structures.py` and/or directly into nf-core/proteinfold.

---

## `update_samplesheet_and_structures.py`

`tools/update_samplesheet_and_structures.py` automates the maintenance loop once you have a MANE missing set. It removes entries already satisfied by the AlphaFold canonical download, reuses those canonical structures whenever possible, and can be run before nf-core/proteinfold (to harvest/prune) and after nf-core completes (to ingest predictions and refresh the remaining missing set).

### Configuration inputs
Non-runtime paths still live in `config.yaml`:

| Field | Description |
| --- | --- |
| `paths.<env>.predicted_dir` | Absolute path to the nf-core/proteinfold run output (predicted ENSP PDBs from nf-core). |
| `paths.<env>.predicted_relpath` | Folder (relative to `--samplesheet-folder`) where the copied predicted ENSP bundle lives (default `predicted/`). |
| `paths.<env>.missing_relpath` | Folder (relative) that will hold `missing/samplesheet.csv` and `missing/fasta/` for the next nf-core proteinfold run (default `missing/`). |
| `paths.<env>.retrieved_relpath` | Folder (relative) where canonical structures reused from the AlphaFold DB canonical PDBs are stored (default `retrieved/`). |
| `paths.<env>.final_bundle_relpath` | Folder (relative) containing the merged `pdbs/` + `samplesheet.csv` to be passed to `oncodrive3d build-datasets` (default `final_bundle/`). |
| `paths.<env>.mane_missing_path` | Absolute path to `mane_missing.csv` (UniProt ↔ RefSeq mapping used for canonical reuse). Produced when running `oncodrive3d build-datasets --mane` or default (see prerequisites). |
| `paths.<env>.mane_summary_path` | Absolute path to the MANE summary file used to map ENSP ↔ gene symbols. Downloaded automatically during `oncodrive3d build-datasets --mane` or default (see prerequisites). |
| `paths.<env>.cgc_list_path` | Absolute path to the Cancer Gene Census TSV (optional; only needed for CGC prioritisation). Download available from the CGC website (registration required). |

Runtime arguments override or complement the config file; `config.yaml` no longer stores `samplesheet_folder` or runtime flags.

### Usage

```bash
python -m tools.preprocessing.update_samplesheet_and_structures \
    --samplesheet-folder /path/to/mane_missing/data \
    [--predicted-dir   /path/to/nfcore/pdbs] \
    [--canonical-dir   /path/to/af_canonical_pdbs]
```

Arguments:

- `--samplesheet-folder` (**required**): Folder created by `prepare_samplesheet.py` (contains `samplesheet.csv` + `fasta/`).
- `--config-path`: YAML with path templates for your local/cluster environments (default `config.yaml`).
- `--predicted-dir`: Directory containing nf-core/proteinfold PDBs to ingest; omit it if you only want to reuse AF canonical structures in this pass.
- `--canonical-dir`: Directory with the AlphaFold canonical PDBs (typically `<build_folder>/pdb_structures` from `oncodrive3d build-datasets --mane` or the default build). Required to reuse canonical structures; `--mane_only` builds alone do not download these files.
- `--max-workers`: Parallel workers for canonical indexing (default = all cores).
- `--enable-canonical-reuse/--disable-canonical-reuse`: Toggle canonical harvesting (enabled by default when `--canonical-dir` is provided).
- `--filter-long-sequences/--no-filter-long-sequences`: Whether to drop long proteins from the nf-core input (default enabled).
- `--max-sequence-length`: Length cutoff applied when filtering (default `2700` residues).

> [!NOTE]
> Provide at least one of `--predicted-dir` or `--canonical-dir`; otherwise the script has nothing to consume.

### Workflow

1. **Sync nf-core predictions (optional)** – copy the raw nf-core/proteinfold PDBs into `predicted/` and emit an accompanying `samplesheet.csv`.
2. **Refresh the missing set** – copy FASTAs for still-missing ENSPs into `missing/fasta/`, rebuild `missing/samplesheet.csv`, and annotate with CGC metadata plus sequence lengths.
3. **Reuse canonical structures (optional)** – if `--canonical-dir` is supplied, index the AlphaFold canonical bundle, map ENSPs via RefSeq IDs, and copy matching PDBs into `retrieved/`, pruning them from the `missing/` tree.
4. **Filter long sequences** – optionally drop proteins exceeding `--max-sequence-length`, logging them in `missing/samplesheet_removed_long.csv`.
5. **Merge bundles** – combine `predicted/` and `retrieved/` into `final_bundle/`, the directory you will pass to `oncodrive3d build-datasets --custom_mane_*`.

> [!NOTE]
> You can run the script twice per iteration: first with only `--canonical-dir` (before nf-core predictions are available) and later with `--predicted-dir` once the nf-core run finishes to integrate the new PDBs.

### Outputs

After each run, `<samplesheet_folder>` contains:

- `predicted/`
  - `samplesheet.csv`
  - `pdbs/ENSP....alphafold.pdb` (copied nf-core predictions)
- `missing/`
  - `samplesheet.csv` (prioritised/filtered for the next nf-core submission)
  - `samplesheet_removed_long.csv` (only when filtering drops entries)
  - `fasta/ENSP....fasta`
- `retrieved/`
  - `samplesheet.csv`
  - `pdbs/...` (canonical structures reused from the AlphaFold DB canonical set)
- `final_bundle/`
  - `samplesheet.csv` (deduplicated union of predicted + retrieved ENSPs; ready for `oncodrive3d build-datasets --custom_mane_metadata_path`)
  - `pdbs/...` (ready for `oncodrive3d build-datasets --custom_mane_pdb_dir`)

---

## End-to-end loop

0. **Bootstrap the datasets** – run `oncodrive3d build-datasets --mane_only` to generate the MANE-only baseline and mapping files, and run `oncodrive3d build-datasets --mane` (or a default `build-datasets`) if you also wish to retrieve structures from the canonical AlphaFold download (recommended).
1. **Prepare the missing set** – run `prepare_samplesheet.py`.
2. **Harvest canonical matches (optional)** – invoke `update_samplesheet_and_structures.py` with `--canonical-dir` to reuse any AlphaFold canonical structures and shrink the missing set before prediction.
3. **Predict the remaining structures** – run nf-core/proteinfold on `missing/samplesheet.csv` + `missing/fasta/`.
4. **Ingest predictions** – re-run `update_samplesheet_and_structures.py` with `--predicted-dir` (and optionally `--canonical-dir` again) to fold new PDBs into `predicted/` and refresh the missing set.
5. **Rebuild Oncodrive3D datasets** – point `oncodrive3d build-datasets --mane_only` at `final_bundle/pdbs` (`--custom_mane_pdb_dir`) and `final_bundle/samplesheet.csv` (`--custom_mane_metadata_path`) so every subsequent `oncodrive3d run` benefits from the extended MANE coverage.
