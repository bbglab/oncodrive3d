# Annotation Build & Plotting Workflow

Oncodrive3D ships with a plotting pipeline that turns the clustering results into visualizations and annotated tables. They help you interpret **why** a gene or residue scored a significant 3D-clustering signal, judge whether that signal looks biologically plausible or artifactual, and produce downstream assets (summary panels, per-gene tracks, association charts, annotated CSVs) for follow-up.

Annotations are built once per dataset via `oncodrive3d build-annotations`, then reused by `oncodrive3d plot` for any cohort.

---

## Prerequisites

Building the annotation bundle (`oncodrive3d build-annotations`) is the demanding step and needs all of the following. Plotting (`oncodrive3d plot`) needs only items 1 and 2, plus the annotation bundle and the outputs of an `oncodrive3d run` (listed under [Generating Plots](#generating-plots)).

1. **Datasets** – `oncodrive3d build-datasets` (or `build-datasets --mane_only`) must have been run already. Both steps read from this folder (`build-annotations` uses `pdb_structures` and `seq_for_mut_prob.tsv`; `plot` uses `confidence.tsv` and `seq_for_mut_prob.tsv`).
2. **Python environment** – use the same environment (e.g., `uv` virtualenv) that you rely on for the CLI.
3. **PDB_Tool binary** – `oncodrive3d build-annotations` invokes the [PDB_Tool](https://github.com/realbigws/PDB_Tool) executable named `PDB_Tool` on `$PATH` to compute per-residue solvent accessibility and secondary structure. See **Installing PDB_Tool** below for a recipe.
4. **Internet access** – required to download Pfam annotations, UniProt features, and (unless `--ddg_dir` is set) RaSP ΔΔG predictions.
5. **Disk space** – annotation folders contain many files; keep several GB free.

### Installing PDB_Tool

[PDB_Tool](https://github.com/realbigws/PDB_Tool) compiles from source. If your conda env has no C/C++ toolchain, install one first:

```bash
conda install -c conda-forge gxx make
```

Then build and put it on `$PATH`:

```bash
git clone https://github.com/realbigws/PDB_Tool.git
cd PDB_Tool
make -C source_code
# The Makefile drops the binary at the repo root. Symlink into the active conda env so it resolves via `which PDB_Tool`:
ln -s "$(pwd)/PDB_Tool" "$CONDA_PREFIX/bin/PDB_Tool"
```

---

## Building the Annotation Bundle

Run once per dataset (or whenever you update AlphaFold structures):

```bash
# Default (Homo sapiens, public RaSP ΔΔG):
oncodrive3d build-annotations -d <build_folder> -o <annot_folder>

# With a MANE-built dataset:
oncodrive3d build-annotations -d <mane_build_folder> -o <annot_folder>

# Mouse with custom ΔΔG predictions (omit --ddg_dir to skip the ΔΔG step):
oncodrive3d build-annotations -d <build_folder> -o <annot_folder> -s mouse --ddg_dir <ddg_path>
```

See `oncodrive3d build-annotations --help` for all options.

**Worth knowing:**

- ΔΔG predictions default to the public RaSP bundle (computed against the canonical AlphaFold v4 human proteome, not the MANE bundle). For datasets built with a different AF version, residue-level mismatches are filtered during validation (see `--ddg_mismatch_threshold` below).
- `--ddg_dir` overrides the default download with a folder of RaSP-style CSVs (columns `variant` and `score_ml`; UniProt accession auto-detected anywhere in the filename, any separator). For non-human organisms the public bundle doesn't apply, so without `--ddg_dir` the ΔΔG step is skipped with a warning; other annotations build normally. To generate predictions yourself on CPUs, see [bbglab/rasp_cpu](https://github.com/bbglab/rasp_cpu).
- `--ddg_mismatch_threshold` (default `0.1`) drops a protein if its wild-type residues disagree with the canonical UniProt sequence above this fraction. Set to `1.0` to disable the WT-mismatch check (positions outside the canonical sequence still drop the protein).
- If `--output_dir` exists and isn't empty, you're prompted before its contents are cleaned (excluding `log/`); pass `--yes` to auto-confirm.

The command assembles three annotation tracks:

- **Stability change (ΔΔG)** – RaSP predictions are downloaded (human) or read from `--ddg_dir`, then validated against the canonical sequence from `seq_for_mut_prob.tsv`; proteins failing validation are dropped with a warning. Positions with no prediction are kept as `NaN` (not `0.0`) so plots show gaps and the logistic regression uses only real measurements.
- **PDB features** – AlphaFold structures are run through `PDB_Tool` into `pdb_tool_df.tsv`, with residue-level secondary structure (`SSE`) and relative accessibility (`pACC`).
- **Domains and sites** – Pfam coordinates (Ensembl BioMart) go to `pfam.tsv`; UniProt DOMAIN/PTM/SITE/MOTIF/MEMBRANE features (EMBL-EBI Proteins API) go to `uniprot_feat.tsv`.

### Output Layout

After a successful run the annotation folder contains:

```text
annotations/
├── pdb_tool_df.tsv
├── pfam.tsv
├── uniprot_feat.tsv
├── stability_change/        # optional; absent for mouse builds without --ddg_dir
│   └── <UNIPROT>_ddg.json
└── log/
```

Keep this directory around: `oncodrive3d plot` reads the tables above and merges in ΔΔG values when `stability_change/` is present, otherwise the ΔΔG track is omitted from per-gene plots and association analyses.

---

## Generating Plots

Once you have:

- Gene-level results (`<cohort>.3d_clustering_genes.csv`),
- Residue-level results (`<cohort>.3d_clustering_pos.csv`),
- Processed mutations (`<cohort>.mutations.processed.tsv`),
- Missense probability dictionary (`<cohort>.miss_prob.processed.json`),
- Processed sequence dataframe (`<cohort>.seq_df.processed.tsv`),
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
  --cohort COHORT
```

See `oncodrive3d plot --help` for all options.

**Worth knowing:**

- `--maf_path` is the **processed missense-only** TSV (`<cohort>.mutations.processed.tsv`) from `oncodrive3d run`. `--maf_for_nonmiss_path` is optional and takes the **original** MAF (before processing); supply it to enable the non-missense track. All other input files (gene/pos results, `--miss_prob_path`, `--seq_df_path`, `--datasets_dir`, `--annotations_dir`) must come from that same `oncodrive3d run` invocation; mismatch yields empty plots or missing-track errors.
- `--lst_summary_tracks` / `--lst_gene_tracks` accept comma-separated track names; pair them with `--lst_*_hratios` to redistribute vertical space.

The command produces:

- A **summary plot** of per-gene mutation counts, cluster residues, and score distributions.
- **Per-gene plots** overlaying the requested tracks (observed vs expected mutation counts, missense probabilities, clustering scores, PAE/pLDDT, ΔΔG, Pfam/UniProt domains, PTMs, membrane regions, motifs). Tracks not available for a gene are dropped automatically.
- **Annotated tables** (`<cohort>.3d_clustering_pos.annotated.csv` and `<cohort>.uniprot_feat.tsv`), merging the positional results with disorder (pLDDT), PDB features, transcript metadata, and UniProt domains. Written only when `--output_csv` is passed.
- **Association plots** (optional); see the "Association Analyses" section below.

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

For interactive-ready 3D views, the separate `oncodrive3d chimerax-plot` command renders PNG snapshots (plus `.defattr` attribute files) under `<output_dir>/<cohort>.chimerax/`. It reuses the gene/position CSVs from `oncodrive3d run`, the datasets directory (for AlphaFold structures), and the processed sequence dataframe. Each snapshot colours the AlphaFold model by a mutation or clustering metric (mutations in residue, mutations in volume, clustering score, log clustering score), highlighting the mutated or cluster residues as spheres. The Nextflow pipeline exposes the same functionality through the `chimerax_plot` flag.

> [!NOTE]
> **ChimeraX must be installed separately.** The framework was tested with **ChimeraX 1.6.1** (`ucsf-chimerax_1.6.1ubuntu20.04_amd64.deb` from [UCSF older releases](https://www.cgl.ucsf.edu/chimerax/older_releases.html); newer releases should also work).
>
> If instead you run Oncodrive3D through the provided `chimerax` or `full` Docker image, ChimeraX is already included, so no separate install is needed.
>
> By default the command looks for the executable at `/usr/bin/chimerax`; pass `--chimerax_bin` if yours is elsewhere.

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

- `--pixel_size` controls resolution: smaller values produce larger images (default `0.08`).
- `--cluster_ext` displays extended clusters (mutations that contribute to but don't directly form significant clusters).
- `--af_version` is auto-detected from the structures in the datasets directory, so you normally don't set it. It's used only as a tiebreaker when the dataset contains more than one AlphaFold version, or as a fallback if none is detected (default `6`).
- `--spheres` / `--no-spheres` (default on) highlights residues as spheres: mutated residues on the base plots, cluster residues on the `*_clusters` plots. With `--no-spheres` the base plots are cartoon-only while the `*_clusters` plots still mark the clusters.
- `--cluster_markers` (default off) adds translucent volume bubbles on the cluster residues in the `*_clusters` plots.
- `--non_mutated_color` (default `gray`) and `--text_color` (default `black`) set the colour of the non-mutated cartoon and of the title / color-bar label (any ChimeraX colour name or hex).
- `--transparent_bg` / `--no-transparent_bg` (default on) saves images with a transparent background; pass `--no-transparent_bg` for a white background.

---

## Troubleshooting & Tips

- **PDB_Tool missing** – install the binary or adjust `$PATH`. The build step logs the exact command being executed, making it easier to diagnose permission issues.
- **Annotation mismatches** – plots rely on UniProt IDs matching between the run outputs and the annotation tables. Make sure you pass the same datasets directory used during `build-datasets`.
- **Association plots without data** – if a gene lacks both clustered and non-clustered residues after filtering, it is skipped from the logistic regression. The log file will note “There aren’t any relationship to plot”.
