# Annotation Build & Plotting Workflow

Oncodrive3D ships with a plotting pipeline that turns the clustering results into analysis-ready visualizations and enriched tables. Its goal is threefold: 
1. help you interpret **why** a gene or residue achieved a significant 3D-clustering signal by overlaying structural/functional context, 
2. highlight diagnostic cues that reveal whether a signal looks biologically plausible or could be an artifact,
3. provide downstream-analysis assets such as summary panels, per-gene tracks, volcano/log-odds association charts, and annotated CSVs to guide follow-up experiments or reporting. 

To keep the runtime reasonable, structural and functional annotations are generated once per dataset via `oncodrive3d build-annotations`, then reused by `oncodrive3d plot` for any cohort.

This document describes the prerequisites, the intermediate files that are produced, and how to customize the plots.

---

## Prerequisites

1. **Datasets** – `oncodrive3d build-datasets` (or `build-datasets --mane_only`) must have been run already; the plotting stage reads `datasets/pdb_structures`, `seq_for_mut_prob.tsv`, and `confidence.tsv`.
2. **Python environment** – use the same environment (e.g., `uv` virtualenv) that you rely on for the CLI.
3. **PDB_Tool binary** – `oncodrive3d build-annotations` invokes the [PDB_Tool](https://github.com/realbigws/PDB_Tool) executable named `PDB_Tool` on `$PATH` to compute per-residue solvent accessibility and secondary structure. Install it or wrap it in a container accessible from the command line.
4. **Internet access** – required to download Pfam annotations, UniProt features, and (if `--ddg_dir` is not set) RaSP ΔΔG predictions. You can point `--ddg_dir` to precomputed RaSP files to skip the download step or to provide in-house predicted scores.
5. **Disk space** – annotation folders contain many files; keep several GB free.

---

## Building the Annotation Bundle

Run once per dataset (or whenever you update AlphaFold structures):

```bash
# Default (Homo sapiens, public RaSP ΔΔG):
oncodrive3d build-annotations -d <build_folder> -o <annot_folder>

# With a MANE-built dataset:
oncodrive3d build-annotations -d <mane_build_folder> -o <annot_folder>

# Mouse (custom ΔΔG predictions required):
oncodrive3d build-annotations -d <build_folder> -o <annot_folder> -s mouse --ddg_dir <ddg_path>
```

See `oncodrive3d build-annotations --help` for all options.

**Worth knowing:**

- `--ddg_dir` must point to a folder of RaSP-style CSVs (columns `variant` and `score_ml`; UniProt accession auto-detected anywhere in the filename, any separator). Required for mouse if ΔΔG tracks are desired — the public RaSP bundle is human-only; omitting it for mouse skips ΔΔG with a warning.
- `--ddg_mismatch_threshold` (default `0.1`) drops a protein if its wild-type residues disagree with the canonical UniProt sequence above this fraction. Set to `1.0` to disable the WT-mismatch check (positions outside the canonical sequence still drop the protein).
- `--output_dir` is cleaned at each run unless `--yes` is passed.

What happens internally (`scripts/plotting/build_annotations.py` and helpers):

1. **Cleanup** – the target directory is emptied (except for `log/`) unless `--yes` is provided.
2. **Stability change (ΔΔG)** – RaSP predictions are downloaded (human) or read from `--ddg_dir`. Each protein is parsed into `{position: {ALT: ddg}}` (averaging across fragments) and validated against the canonical sequence from `seq_for_mut_prob.tsv`; proteins failing validation are dropped with a warning. Within a kept protein, positions with no prediction surface as `NaN` (not `0.0`) so plots show gaps, the annotated CSV distinguishes "no data" from "neutral mutation", and the logistic regression restricts itself to real measurements.
3. **PDB features** – AlphaFold structures are decompressed and sent through `PDB_Tool`, producing `.feature` files that are then parsed into `pdb_tool_df.tsv` with residue-level secondary structure (`SSE`) and relative accessibility (`pACC`).
4. **Pfam domains** – Pfam coordinates are pulled from the Ensembl BioMart archive plus the Pfam ID database, merged with Oncodrive3D’s sequence metadata, and written to `pfam.tsv`.
5. **UniProt features** – the EMBL-EBI Proteins API supplies DOMAIN/PTM/SITE/MOTIF/MEMBRANE annotations. They are normalized, merged with Pfam entries, and stored in `uniprot_feat.tsv`.

### Output Layout

After a successful run the annotation folder contains:

```
annotations/
├── pdb_tool_df.tsv
├── pfam.tsv
├── uniprot_feat.tsv
├── stability_change/        # optional; absent for mouse builds without --ddg_dir
│   └── <UNIPROT>_ddg.json
└── log/
```

Keep this directory around — `oncodrive3d plot` reads the tables above and merges in ΔΔG values when `stability_change/` is present, otherwise the ΔΔG track is omitted from per-gene plots and association analyses.

---

## Generating Plots

Once you have:

- Gene-level results (`<cohort>.3d_clustering_genes.csv`),
- Residue-level results (`<cohort>.3d_clustering_pos.csv`),
- Processed mutations (`<cohort>.mutations.processed.tsv`),
- Missense probability dictionary (`<cohort>.miss_prob.processed.json`),
- `seq_for_mut_prob.tsv`,
- Built datasets (`datasets/`) and annotations (`annotations/`),

call:

```bash
oncodrive3d plot \
  --gene_result_path output/COHORT/COHORT.3d_clustering_genes.csv \
  --pos_result_path output/COHORT/COHORT.3d_clustering_pos.csv \
  --maf_path output/COHORT/COHORT.mutations.processed.tsv \
  --miss_prob_path output/COHORT/COHORT.miss_prob.processed.json \
  --seq_df_path output/COHORT/COHORT.seq_df.processed.tsv \
  --datasets_dir /path/to/datasets \
  --annotations_dir /path/to/annotations \
  --output_dir plots/COHORT \
  --cohort COHORT \
  --maf_for_nonmiss_path original_input.maf \
  --lst_gene_tracks miss_count,miss_prob,score,clusters,ddg,disorder,pacc,ptm,site,sse,pfam
```

See `oncodrive3d plot --help` for all options.

**Worth knowing:**

- All input files (gene/pos results, `--maf_path`, `--miss_prob_path`, `--seq_df_path`) must come from the same `oncodrive3d run` invocation; `--datasets_dir` and `--annotations_dir` must be the same pair used during that run. Any mismatch yields empty plots or missing-track errors.
- `--maf_for_nonmiss_path` is optional and points to the **original** MAF (not the processed one); supplying it enables the non-missense track.
- `--lst_summary_tracks` / `--lst_gene_tracks` accept comma-separated track names; pair them with `--lst_*_hratios` to redistribute vertical space.

During execution (`scripts/plotting/plot.py`):

1. Results are filtered to genes requested by the user, and the corresponding entries are sliced out of the sequence/annotation tables.
2. A **summary plot** shows per-gene mutation counts, cluster residues, and score distributions.
3. **Per-gene plots** combine multiple tracks: observed vs expected mutation counts, missense probabilities, clustering scores, PAE/pLDDT, ΔΔG, Pfam/UniProt annotations, PTMs, membrane regions, motifs, etc. Tracks that are not available for a gene are automatically removed.
4. Annotated tables are built by merging the positional results with disorder (pLDDT), PDB features, transcript metadata, and UniProt domains. They are saved as `<cohort>.3d_clustering_pos.annotated.csv` plus `<cohort>.uniprot_feat.tsv`.
5. **Association plots** (optional) – see the “Association Analyses” section below for details on how the logistic-regression statistics are generated and visualized (volcano, per-gene volcano, and log-odds panels).

---

### Association Analyses

The optional association module quantifies how strongly specific annotations track with significant clusters:

1. **Input preparation** – residues with non-zero missense probability inherit standardized predictors: structural metrics (pLDDT, ΔΔG, surface exposure), categorical features (Pfam/UniProt/PTM/motif dummies), and the expected missense probability itself.
2. **Univariate logistic regressions** – for each gene and each predictor, the pipeline fits `logit(C ~ feature)` where `C` is the binary cluster label. This yields log-odds, standard errors, and raw p-values that are stored in `<cohort>.logreg_result.tsv`.
3. **Visualization** – the statistics above feed three plot types under `<cohort>.associations_plots/`:
   - A cohort-wide volcano plot highlighting annotations with the most extreme log-odds and p-values.
   - Per-gene mini volcano plots to inspect feature enrichments gene by gene.
   - Log-odds strip charts with 95% confidence intervals to visualize effect sizes.

Only raw p-values are provided; apply your preferred multiple-testing correction (e.g., BH-FDR) before drawing conclusions about specific features.

---

### ChimeraX 3D Snapshots

For interactive-ready 3D views, Oncodrive3D exposes a separate `oncodrive3d chimerax-plot` command. It takes the same gene/position-level CSVs produced by `oncodrive3d run`, plus the datasets directory (for AlphaFold structures) and the processed sequence dataframe, and renders PNG snapshots together with `.defattr` files under `<output_dir>/<cohort>.chimerax/`. Each snapshot shows the AlphaFold model with significant clusters highlighted, optional extended clusters, pLDDT coloring, and sample counts. Provide the path to your ChimeraX installation via `--chimerax_bin` or rely on the default `/usr/bin/chimerax`. The Nextflow pipeline exposes the same functionality through the `chimerax_plot` flag.

Example:

```bash
oncodrive3d chimerax-plot \
  --gene_result_path output/COHORT/COHORT.3d_clustering_genes.csv \
  --pos_result_path output/COHORT/COHORT.3d_clustering_pos.csv \
  --datasets_dir /path/to/datasets \
  --seq_df_path output/COHORT/COHORT.seq_df.processed.tsv \
  --output_dir plots/COHORT \
  --cohort COHORT \
  --chimerax_bin /opt/ChimeraX/bin/ChimeraX \
  --max_n_genes 20 \
  --pixel_size 0.1 \
  --cluster_ext
```

See `oncodrive3d chimerax-plot --help` for all options.

**Worth knowing:**

- `--chimerax_bin` defaults to `/usr/bin/chimerax`; override it if ChimeraX is installed elsewhere or running in a container.
- `--pixel_size` controls resolution — smaller values produce larger images (default `0.08`).
- `--cluster_ext` displays extended clusters (mutations that contribute to but don't directly form significant clusters).
- `--af_version` must match the AlphaFold build behind your datasets (default `4`) so structure filenames resolve correctly.

---

## Troubleshooting & Tips

- **PDB_Tool missing** – install the binary or adjust `$PATH`. The build step logs the exact command being executed, making it easier to diagnose permission issues.
- **Annotation mismatches** – plots rely on UniProt IDs matching between the run outputs and the annotation tables. Make sure you pass the same datasets directory used during `build-datasets`.
- **Association plots without data** – if a gene lacks both clustered and non-clustered residues after filtering, it is skipped from the logistic regression. The log file will note “There aren’t any relationship to plot”.
