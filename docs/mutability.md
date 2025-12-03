# Mutability-Aware Mutation Probabilities

Oncodrive3D normally derives per-residue missense probabilities from a **mutational profile** (192 trinucleotide contexts) that reflects the cohort-wide substitution spectrum. This assumption breaks when sequencing depth varies substantially across loci—as is typical for ultra-deep Duplex sequencing assays used to profile normal tissues. In those cases you can feed Oncodrive3D a **site-level mutability table** so that the background model reflects the actual opportunity to call mutations in each coding base. The table itself captures both the local sequencing depth/callability and the substitution bias, so no external mutational profile is required.

This document explains the required inputs, how the mutability data are consumed by the pipeline, and how to interpret the additional diagnostics that appear in the outputs.

---

## Input Requirements

### 1. Oncodrive3D sequence datasets

`oncodrive3d build-datasets` must have been executed beforehand. The run produces `seq_for_mut_prob.tsv`, which stores—for every UniProt fragment used during clustering—the reference coding DNA sequence, per-base trinucleotide contexts, exon genomic coordinates, and flags indicating whether those coordinates come from the Proteins API (`Reference_info == 1`).  
Only entries with `Reference_info == 1` can be paired with site-level mutability values, because their exon coordinates are known precisely. Genes that lack this metadata will be skipped when mutability mode is enabled.

### 2. Site-level mutability table

Provide a tab-delimited file describing every callable coding site. The file must be **bgzip-compressed and indexed with tabix**, because Oncodrive3D uses `pysam/tabix` to query genomic windows gene by gene (`scripts/run/mutability.py`). Each line should include:

- Chromosome (with or without a `chr` prefix; configurable),
- Genomic position (1-based),
- Reference nucleotide,
- Alternate nucleotide,
- Mutability score (float between 0 and 1),
- Optionally, an element identifier (e.g., transcript, exon, panel tile) if you need to disambiguate overlapping coordinates.

The caller is responsible for generating these values (e.g., by combining coverage data, callability masks, and the mutational profile). Oncodrive3D simply treats the mutability number as the probability that a given ref→alt substitution could be detected at that site.

### 3. JSON configuration

Because mutability tables vary across cohorts, Oncodrive3D expects a minimal JSON file describing how to parse yours. All offsets are **0-based column indices**:

```json
{
  "file": "/path/to/all_samples.mutability_per_site.tsv.gz",
  "chr_prefix": "chr",
  "chr": 0,
  "pos": 1,
  "ref": 2,
  "alt": 3,
  "mutab": 4,
  "element": 5,
  "extra": null
}
```

Required keys:

| Key | Description |
| --- | --- |
| `file` | Path to the bgzip-compressed, tabix-indexed TSV. |
| `chr_prefix` | Prefix to prepend when querying tabix (`""` if your file stores bare chromosome numbers). |
| `chr` / `pos` | Column indices for chromosome and genomic position. |
| `ref` / `alt` | Column indices for the reference and alternate nucleotide; use `null` if unavailable. |
| `mutab` | Column index containing the mutability score to use. |

Optional keys:

| Key | Description |
| --- | --- |
| `element` | Column storing an element identifier; used to double-check that rows belong to the current transcript. |
| `extra` | Reserved for additional metadata; currently unused but accepted to keep configs forward-compatible. |

---

## Running Oncodrive3D with Mutability Data

1. Ensure the mutability TSV is bgzip-compressed and tabix-indexed:
   ```bash
   bgzip -@ 8 all_samples.mutability_per_site.tsv
   tabix -s 1 -b 2 -e 2 all_samples.mutability_per_site.tsv.gz
   ```
   Adjust the `-s/-b/-e` arguments to match your column order.

2. Create the JSON configuration describing the column offsets (see above).

3. Launch `oncodrive3d run` with the `--mutability_config_path` flag:
   ```bash
   oncodrive3d run \
     -i cohort.maf \
     -d /path/to/datasets \
     -o output/cohort \
     --mutability_config_path configs/mutability.json \
     --thr_mapping_issue 0.2
   ```
   - When this flag is present, the mutational profile (`--mut_profile_path`) is ignored because the site-level table already encodes both sequencing-depth effects and substitution preferences.
   - Genes whose sequences lack `Reference_info == 1` are reported with status `No_mutability`.

Behind the scenes (`scripts/run/clustering.py` and `scripts/run/miss_mut_prob.py`):

- The JSON config simply tells Oncodrive3D how to read your table. Once loaded, tabix queries pull the mutability entries overlapping each gene fragment’s exon coordinates.
- Every coding base of the fragment inherits a **coverage-aware mutation opportunity**: the software respects transcript strand orientation, fills any gaps with zero probability, and builds an internal map of “position → probability for each possible substitution”.
- Those base-level opportunities are collapsed to per-codon / per-residue missense probabilities, which replace the profile-based ones during simulations and scoring.
- The final per-residue vectors are serialized to `<cohort>.miss_prob.processed.json` so that downstream plotting, QC, and association analyses use the exact same background model.

---

## Diagnostics in the Output

When mutability mode is active, additional columns and statuses appear in the gene-level report (`<cohort>.3d_clustering_genes.csv`):

| Field / Status | Meaning |
| --- | --- |
| `Mut_zero_mut_prob` | Number of input mutations that landed on residues with zero mutability in the provided table. |
| `Pos_zero_mut_prob` | Specific residue positions (1-based protein coordinates) whose mutability was zero. |
| `Ratio_mut_zero_prob` | Fraction of the gene’s mutations with zero mutability (used internally to compare against `--thr_mapping_issue`). |
| `Status = Mut_with_zero_prob` | Raised when the ratio above exceeds the selected threshold; the gene is skipped from clustering. |
| `Status = No_mutability` | No mutability data were available for the gene (typically because the transcript lacks high-quality exon coordinates). |

If `--thr_mapping_issue` is set to `1`, genes are never dropped due to zero-miss-probability positions; only the offending mutations are filtered out. Lower thresholds (default `0.1`) provide stricter safeguards for poorly covered regions.
