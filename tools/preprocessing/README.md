# MANE Preprocessing Toolkit

These tools are provided for users that want to run Oncodrive3D using only structures associated to MANE Select transcripts while still covering as many proteins as possible. In fact, AlphaFold database MANE download bundle does not yet contain structures for every MANE Select transcript. Oncodrive3D relies on that bundle when building its datasets, so missing structures translate into genes that cannot be analyzed. The scripts in `tools/preprocessing/` close this gap.  

- `prepare_samplesheet.py` scans the full MANE release and emits `samplesheet.csv` plus per-ENSP FASTAs for every MANE structure that is absent from the AlphaFold MANE download.
- `update_samplesheet_and_structures.py` removes the MANE entries already covered by the AlphaFold canonical bundle (if provided), reuses those canonical structures when available, and copy nf-core predictions into the custom bundle while pruning fulfilled entries from the next proteinfold run input files (the `samplesheet.tsv` and the corresponding `fasta/`).

Together they allow to iteratively update the MANE structures and feed them back into `oncodrive3d build-datasets --custom_mane_pdb_dir <path/to/final_bundle/pdbs> --custom_mane_metadata_path  <path/to/final_bundle/samplesheet.csv>`.
If you want to see the whole process in context, jump to the [Example](#example) under [End-to-end loop](#end-to-end-loop).

## Installation
Requires:
- Python 3.10+
- Oncodrive3D datasets built (a `--mane_only` run for the MANE-only baseline plus, if needed, a standard `--mane` or default run for the canonical bundle). These builds supply `mane_refseq_prot_to_alphafold.csv`, `mane_summary.txt.gz`, `mane_missing.csv`, and (optionally) the canonical AlphaFold PDB bundle.

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
    --mane-dataset-dir /path/to/o3d_mane_only_dataset \
    --output-dir   /path/to/mane_missing/data
```

Key options:

- `--mane-dataset-dir/-d` (**required**): Directory with the MANE-only dataset (`oncodrive3d build-datasets --mane_only`) containing `mane_refseq_prot_to_alphafold.csv` and `mane_summary.txt.gz`.
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

### Purpose
Automates the maintenance loop once you have a MANE missing set. It removes entries already satisfied by the AlphaFold canonical bundle, reuses those canonical structures whenever possible, and can be run before nf-core/proteinfold (to harvest/prune) and after nf-core completes (to ingest predictions and refresh the remaining missing set).

### Configuration inputs
Non-runtime paths still live in `config.yaml`:

| Field | Description |
| --- | --- |
| `paths.samplesheet_relpath` | Template for the samplesheet file (default `{samplesheet_folder}/samplesheet.csv`). |
| `paths.fasta_relpath` | Template for the FASTA directory (default `{samplesheet_folder}/fasta`). |
| `paths.missing_relpath` | Template for the next nf-core batch inputs (default `{samplesheet_folder}/missing`). |
| `paths.predicted_relpath` | Template for the nf-core predictions bundle (default `{samplesheet_folder}/predicted`). |
| `paths.retrieved_relpath` | Template for canonical reuse outputs (default `{samplesheet_folder}/retrieved`). |
| `paths.final_bundle_relpath` | Template for the merged bundle passed to `oncodrive3d build-datasets` (default `{samplesheet_folder}/final_bundle`). |

### Usage

```bash
python -m tools.preprocessing.update_samplesheet_and_structures \
    --samplesheet-folder /path/to/mane_missing/data \
    --mane-dataset-dir /path/to/mane_only_dataset \
    [--canonical-dir   /path/to/af_canonical_pdbs] \
    [--predicted-dir   /path/to/nfcore/pdbs]
```

Arguments:

- `--samplesheet-folder` (**required**): Absolute path to the folder created by `prepare_samplesheet.py` (contains `samplesheet.csv` + `fasta/`).
- `--mane-dataset-dir` (**required**): Path to the MANE-only dataset built via `oncodrive3d build-datasets --mane_only`.
- `--canonical-dir`: Path to the dataset including UniProt canonical structures built via `oncodrive3d build-datasets` or `oncodrive3d build-datasets --mane`.
- `--predicted-dir`: Path to the directory containing nf-core/proteinfold PDBs to ingest; omit it if you only want to reuse AF canonical structures in this pass.
- `--cgc-list-path`: Path to the Cancer Gene Census TSV (optional; only needed for CGC prioritisation). Download available from the CGC website (registration required).
- `--max-workers`: Parallel workers for canonical indexing (default = all cores).
- `--enable-canonical-reuse/--disable-canonical-reuse`: Toggle canonical harvesting (enabled by default when `--canonical-dir` is provided).
- `--filter-long-sequences/--no-filter-long-sequences`: Whether to drop long proteins from the nf-core input (default enabled).
- `--max-sequence-length`: Length cutoff applied when filtering (default `2700` residues).
- `--config-path`: YAML with path templates describing where to place predicted/missing/retrieved bundles relative to `--samplesheet-folder` (default `config.yaml`).

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

0. **Initialize the datasets** – run `oncodrive3d build-datasets --mane_only` to generate the MANE-only baseline and mapping files, and run `oncodrive3d build-datasets --mane` (or a default `build-datasets`) if you also wish to retrieve structures from the canonical AlphaFold download (recommended).
1. **Prepare the missing set** – run `prepare_samplesheet.py`.
2. **Harvest canonical matches (optional)** – invoke `update_samplesheet_and_structures.py` with `--canonical-dir` to reuse any AlphaFold canonical structures and shrink the missing set before prediction.
3. **Predict the remaining structures** – run nf-core/proteinfold on `missing/samplesheet.csv` + `missing/fasta/`.
4. **Ingest predictions** – re-run `update_samplesheet_and_structures.py` with `--predicted-dir` (and optionally `--canonical-dir` again) to fold new PDBs into `predicted/` and refresh the missing set.
5. **Rebuild Oncodrive3D datasets** – point `oncodrive3d build-datasets --mane_only` at `final_bundle/pdbs` (`--custom_mane_pdb_dir`) and `final_bundle/samplesheet.csv` (`--custom_mane_metadata_path`) so every subsequent `oncodrive3d run` benefits from the extended MANE coverage.

### Example

1. **Initialize datasets**
   ```bash
   # Activate the main O3D environment
   source .venv/bin/activate

   # Build O3D datasets
   oncodrive3d build-datasets --mane_only   --output_dir <path/to/o3d_datasets-mane_only-date>
   oncodrive3d build-datasets               --output_dir <path/to/o3d_datasets-date>
   ```

2. **Prepare missing set**
   ```bash
   # Activate tools environment
   cd tools/preprocessing
   source .venv/bin/activate

   # Init the MANE missing structures
   python -m tools.preprocessing.prepare_samplesheet \
     --mane-dataset-dir <path/to/o3d_datasets-mane_only-date> \
     --output-dir       <path/to/mane_missing-date>
   ```

3. **Harvest canonical matches (first iteration, optional but recommended)**
   ```bash
   # Retrieve MANE missing structures overlapping sequences of canonical ones
   python -m tools.preprocessing.update_samplesheet_and_structures \
     --samplesheet-folder <path/to/mane_missing-date> \
     --mane-dataset-dir   <path/to/o3d_datasets-mane_only-date> \
     --canonical-dir      <path/to/o3d_datasets-date> \
     --cgc-list-path      <path/to/cgc_list>            # (optional, use to sort the list)
   ```

4. **Predict remaining structures**
   - Feed `<path/to/mane_missing-date>/missing/{samplesheet.csv,fasta/}` to nf-core/proteinfold.

5. **Ingest predictions + canonical reuse**
   ```bash
   # Merge retrieved + predicted structures into a final_bundle
   python -m tools.preprocessing.update_samplesheet_and_structures \
     --samplesheet-folder <path/to/mane_missing-date> \
     --mane-dataset-dir   <path/to/o3d_datasets-mane_only-date> \
     --canonical-dir      <path/to/o3d_datasets-date> \
     --predicted-dir      <path/to/predicted/pdbs>      # (what nf-core/proteinfold produces)
     --cgc-list-path      <path/to/cgc_list>            # (optional, use to sort the list)
   ```

6. **Rebuild MANE-only datasets with the final bundle**
   ```bash
   # Build a new MANE only datasets providing the added structures in the final bundle
   oncodrive3d build-datasets --mane_only \
     --custom_mane_pdb_dir          <path/to/mane_missing-date>/final_bundle/pdbs \
     --custom_mane_metadata_path    <path/to/mane_missing-date>/final_bundle/samplesheet.csv \
     --output_dir                   <path/to/mane_missing-new_date>
   ```
