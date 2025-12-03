# Oncodrive3D Nextflow Pipeline

This pipeline enables running Oncodrive3D in parallel across multiple cohorts using [Nextflow](https://www.nextflow.io/).

## Requirements

1. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) (version `23.04.3` was used for testing).
2. Install and set up either or both:
   - [Singularity](https://sylabs.io/guides/latest/user-guide/installation.html)  
   - [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)  
      Ensure Oncodrive3D is installed in your Conda environment and update the `params` section of the `nextflow.config` file to point to your Conda installation:

         ```groovy
         params {
            ...
            conda_env = '/path/to/conda/environment/with/oncodrive3d' 
            ...
         }
         ```

      Replace `/path/to/conda/environment/with/oncodrive3d` with the path to your Conda environment. Alternatively, you can provide it as a command-line argument.


## Test Run

Run a test to ensure that everything is set up correctly and functioning as expected:

```
nextflow run main.nf -profile test,container --data_dir <build_folder>
```

Replace `<build_folder>` with the path to the Oncodrive3D datasets built in the [building datasets](../README.md#building-datasets) step.
If you prefer to use Conda, replace `container` in the `-profile` argument with `conda`.

## Usage

---

> [!WARNING]
> When using the Nextflow script, ensure that your input files are organized in the following directory structure (you only need either the `maf/` or `vep/` directory):
> 
> ```plaintext
> input/
>   ├── maf/
>   │   └── <cohort>.in.maf
>   ├── vep/
>   │   └── <cohort>.vep.tsv.gz
>   └── mut_profile/
>       └── <cohort>.sig.json
> ```
> 
> - `maf/`: Contains mutation files with the `.in.maf` extension.
> - `vep/`: Contains VEP annotation files with the `.vep.tsv.gz` extension, which include annotated mutations with all possible transcripts.
> - `mut_profile/`: Contains mutational profile files with the `.sig.json` extension.

---

```
Usage: nextflow run main.nf [OPTIONS]

Example of run using VEP output as input and MANE Select transcripts:
  nextflow run main.nf -profile container --data_dir <build_folder> --indir <input> \
                       --vep_input true --mane true
  
Options:
  --indir PATH                    Path to the input directory including the subdirectories 
                                  `maf` or `vep` and `mut_profile`. 
  --outdir PATH                   Path to the output directory. 
                                  Default: run_<timestamp>/
  --cohort_pattern STR            Pattern expression to filter specific files within the 
                                  input directory (e.g., 'TCGA*' select only TCGA cohorts). 
                                  Default: *
  --data_dir PATH                 Path to the Oncodrive3D datasets directory, which includes 
                                  the files compiled during the building datasets step.
                                  Default: ${baseDir}/datasets/
  --max_running INT               Maximum number of cohorts to process in parallel.
                                  Default: 5
  --cores INT                     Number of CPU cores used to process each cohort. 
                                  Default: 10
  --memory STR                    Amount of memory allocated for processing each cohort. 
                                  Default: 70GB
  --vep_input BOLEAN              Use `vep/` subdir as input and select transcripts matching 
                                  the Ensembl transcript IDs in Oncodrive3D built datasets. 
                                  Default: false
  --mane                          Prioritize structures corresponding to MANE transcrips if 
                                  multiple structures are associated to the same gene.
                                  Default: false
  --seed INT                      Seed value for reproducibility.
                                  Default: 128
  --plot BOOL                     Generate summary/gene plots (requires annotations_dir).
                                  Default: false
  --chimerax_plot BOOL            Generate ChimeraX snapshots.
                                  Default: false
  --run_extra_args STR            Extra CLI args forwarded to `oncodrive3d run`
                                  (e.g., "--mutability_config_path /path/to/mut.json").
                                  Default: ""
```

### Optional modules

- `--plot true` – runs the same plotting workflow described in [docs/annotations_plotting.md](../docs/annotations_plotting.md). Requires a pre-built annotations folder (`--annotations_dir`) generated via `oncodrive3d build-annotations`. The pipeline publishes summary plots, per-gene panels, annotated CSVs, UniProt feature tables, and association plots under `${outdir}/${outsubdir}/`.
- `--chimerax_plot true` – submits the optional `oncodrive3d chimerax-plot` step for each processed cohort (see [docs/annotations_plotting.md#chimerax-3d-snapshots](../docs/annotations_plotting.md#chimerax-3d-snapshots)). Ensure the `chimerax` binary is accessible inside the container/Conda env or set `process.ext.args` for `O3D_CHIMERAX_PLOT` (in `nextflow.config`) to include `--chimerax_bin /path/to/ChimeraX` if you need a custom location.

### Mutability-aware runs

The pipeline expects a mutational profile JSON (`mut_profile/<cohort>.sig.json`) for every cohort. To supply an additional mutability configuration (e.g., for heterogeneous-depth cohorts), pass the CLI flag through `--run_extra_args`, for example:

```
nextflow run main.nf \
  --run_extra_args "--mutability_config_path /path/to/mutability.json --thr_mapping_issue 0.2"
```

Oncodrive3D prioritizes `--mutability_config_path` over the mutational profile, so the provided `.sig.json` file can be a generic profile if one is not available; it will be ignored once mutability is enabled. Refer to [docs/mutability.md](../docs/mutability.md) for details on preparing the config and TSV.