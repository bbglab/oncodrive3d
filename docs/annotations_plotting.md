# Annotation Build & Plotting Workflow

Oncodrive3D ships with a plotting pipeline that turns the clustering results into analysis-ready visualizations and enriched tables. Its goal is threefold: (1) help you interpret **why** a gene or residue achieved a significant 3D-clustering signal by overlaying structural/functional context, (2) highlight diagnostic cues that reveal whether a signal looks biologically plausible or could be an artifact, and (3) provide downstream-analysis assets such as summary panels, per-gene tracks, volcano/log-odds association charts, and annotated CSVs to guide follow-up experiments or reporting. To keep the runtime reasonable, structural and functional annotations are generated once per dataset via `oncodrive3d build-annotations`, then reused by `oncodrive3d plot` for any cohort.

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
oncodrive3d build-annotations \
  --data_dir /path/to/datasets \
  --output_dir /path/to/annotations
```

What happens internally (`scripts/plotting/build_annotations.py` and helpers):

1. **Cleanup** – the target directory is emptied (except for `log/`) unless `--yes` is provided.
2. **Stability change (ΔΔG)** – for human datasets the command downloads RaSP predictions (`scripts/plotting/stability_change.py`). There is no public RaSP bundle for mouse, so mouse runs must supply their own predictions with `--ddg_dir` if ΔΔG tracks are desired. In all cases, whichever folder you point to is parsed into JSON maps `{position: {ALT: ddg}}`, averaging overlapping fragments.
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
├── stability_change/
│   └── <UNIPROT>_ddg.json
└── log/
```

Keep this directory around—`oncodrive3d plot` expects all of these files.

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

Key options:

- `--lst_summary_tracks` / `--lst_gene_tracks` – comma-separated track names; defaults cover the most common layouts. Matching `--lst_*_hratios` let you tweak relative heights.
- `--maf_for_nonmiss_path` – optional MAF including non-missense mutations to populate the “nonmiss_count” track.
- `--genes` – restrict plotting to specific gene symbols.
- `--c_genes_only` – plot only significant genes (`C_gene == 1`) unless disabled.
- `--volcano_top_n`, `--summary_fsize_*`, `--gene_fsize_*` – control aesthetics of the volcano/log-odds plots and per-gene panels.
- `--output_csv` / `--output_all_pos` – export annotated CSVs, optionally including non-mutated residues.

During execution (`scripts/plotting/plot.py`):

1. Results are filtered to genes requested by the user, and the corresponding entries are sliced out of the sequence/annotation tables.
2. A **summary plot** shows per-gene mutation counts, cluster residues, and score distributions.
3. **Per-gene plots** combine multiple tracks: observed vs expected mutation counts, missense probabilities, clustering scores, PAE/pLDDT, ΔΔG, Pfam/UniProt annotations, PTMs, membrane regions, motifs, etc. Tracks that are not available for a gene are automatically removed.
4. Annotated tables are built by merging the positional results with disorder (pLDDT), PDB features, transcript metadata, and UniProt domains. They are saved as `<cohort>.3d_clustering_pos.annotated.csv` plus `<cohort>.uniprot_feat.tsv`.
5. **Association plots** (Optional) – see the “Association Analyses” section below for details on how the logistic-regression statistics are generated and visualized (volcano, per-gene volcano, and log-odds panels).

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

## Troubleshooting & Tips

- **PDB_Tool missing** – install the binary or adjust `$PATH`. The build step logs the exact command being executed, making it easier to diagnose permission issues.
- **Annotation mismatches** – plots rely on UniProt IDs matching between the run outputs and the annotation tables. Make sure you pass the same datasets directory used during `build-datasets`.
- **Association plots without data** – if a gene lacks both clustered and non-clustered residues after filtering, it is skipped from the logistic regression. The log file will note “There aren’t any relationship to plot”.
