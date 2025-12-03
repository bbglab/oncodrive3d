# Oncodrive3D

**Oncodrive3D** is a fast and accurate computational method designed to analyze patterns of somatic mutation across tumors, with the goal of identifying **three-dimensional (3D) clusters** of missense mutations and detecting genes under **positive selection**. 

The method leverages **AlphaFold 2-predicted protein structures** and Predicted Aligned Error (PAE) to define residue contacts within the protein's 3D space. When available, it integrates **mutational profiles** to build an accurate background model of neutral mutagenesis. By applying a novel **rank-based statistical approach**, Oncodrive3D scores potential 3D clusters and computes empirical p-values.

[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![docker](https://img.shields.io/docker/v/bbglab/oncodrive3d?logo=docker)](https://hub.docker.com/r/bbglab/oncodrive3d)
[![PyPI - Version](https://img.shields.io/pypi/v/oncodrive3d?logo=pypi)](https://pypi.org/project/Oncodrive3D/)

![Graphical abstract of Oncodrive3D](docs/images/graphical_abstract.png "Oncodrive3D")

---

## Requirements

Before you begin, ensure **Python 3.10 or later** is installed on your system.  
Additionally, you may need to install additional development tools. Depending on your environment, you can choose one of the following methods:

- If you have sudo privileges:

   ```bash
   sudo apt install built-essential
   ```

- For HPC cluster environment, it is recommended to use [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) (or [Mamba](https://mamba.readthedocs.io/en/latest/)):

   ```bash
   conda create -n o3d python=3.10.0
   conda activate o3d
   conda install -c conda-forge gxx gcc libxcrypt clang zlib
   ```


## Installation

- Install via PyPI:

   ```bash
   pip install oncodrive3d
   ```

- Alternatively, you can obtain the latest code from the repository and install it for development with pip:

   ```bash
   git clone https://github.com/bbglab/oncodrive3d.git
   cd oncodrive3d
   pip install -e .
   oncodrive3d --help
   ```

- Or you can use a modern build tool like [uv](https://github.com/astral-sh/uv):

   ```bash
   git clone https://github.com/bbglab/oncodrive3d.git
   cd oncodrive3d
   uv run oncodrive3d --help
   ```

## Building Datasets

This step build the datasets necessary for Oncodrive3D to run the 3D clustering analysis. It is required once after installation or whenever you need to generate datasets for a different organism or apply a specific threshold to define amino acid contacts.

> [!WARNING]
> This step is highly time- and resource-intensive, requiring a significant amount of free disk space and computational power. It will download and process a large amount of data. Ensure sufficient resources are available before proceeding, as insufficient capacity may result in extended runtimes or processing failures.
>
> Reliable internet access is required because AlphaFold structures, Ensembl annotations, Pfam files, and other resources are downloaded on demand during the build.

> [!NOTE]
> The first time that you run Oncodrive3D building dataset step with a given reference genome, it will download it from our servers. By default the downloaded datasets go to`~/.bgdata`. If you want to move these datasets to another folder you have to define the system environment variable `BGDATA_LOCAL` with an export command.

> [!NOTE]
> Human datasets built with the default settings pull canonical transcript metadata from the January 2024 Ensembl archive (release 111 / GENCODE v45). For maximum compatibility, annotate your input variants with the same Ensembl/Gencode release or supply the unfiltered VEP output together with `--o3d_transcripts --use_input_symbols`.

```
Usage: oncodrive3d build-datasets [OPTIONS]

Examples:
  Basic build:
    oncodrive3d build-datasets -o <build_folder>
  
  Build with MANE Select transcripts:
    oncodrive3d build-datasets -o <build_folder> --mane

Options:
  -o, --output_dir PATH           Path to the directory where the output files will be saved. 
                                  Default: ./datasets/
  -s, --organism PATH             Specifies the organism (`human` or `mouse`). 
                                  Default: human
  -m, --mane                      Use structures predicted from MANE Select transcripts 
                                  (applicable to Homo sapiens only).
  -M, --mane_only                 Use only structures predicted from MANE Select transcripts
                                  (applicable to Homo sapiens only).
  -C, --custom_mane_pdb_dir PATH  Path to directory containing custom MANE PDB structures.
                                  Default: None
  -f, --custom_mane_metadata_path Path to a dataframe (typically a samplesheet.csv) including 
                                  Ensembl IDs and sequences of the custom pdbs.
                                  Default: None
  -d, --distance_threshold INT    Distance threshold (Å) for defining residues contacts. 
                                  Default: 10
  -c, --cores INT                 Number of CPU cores for computation. 
                                  Default: All available CPU cores
  --af_version INT                Version of the AlphaFold Protein Structure Database release.
                                  Default: 4
  -y, --yes                       Run without interactive prompts.
  -v, --verbose                   Enables verbose output.
  -h, --help                      Show this message and exit.  
```

For more information on the output of this step, please refer to the [Building Datasets Output Documentation](https://github.com/bbglab/oncodrive3d/tree/master/docs/build_output.md).

> [!TIP]
> ### Increasing MANE Structural Coverage
> To maximize structural coverage of **MANE Select transcripts**, you can [predict missing structures locally and integrate them into Oncodrive3D](tools/preprocessing/README.md) using:
>
> - `tools/preprocessing/prepare_samplesheet.py`: a standalone utility that:
>   - Retrieve the full MANE entries from NCBI.
>   - Identifies proteins missing from the AlphaFold MANE dataset.
>   - Generates:
>     - A `samplesheet.csv` with Ensembl protein IDs, FASTA paths, and optional sequences.
>     - Individual FASTA files for each missing protein.
>
> - `tools/preprocessing/update_samplesheet_and_structures.py`: takes the samplesheet folder produced above and:
>   - Reuses canonical AlphaFold structures when possible to shrink the nf-core input.
>   - Ingests nf-core/proteinfold predictions and keeps the `missing/` set up to date.
>   - Generates a `final_bundle/samplesheet.csv` plus `final_bundle/pdbs/`, ready to be passed to `oncodrive3d build-datasets --mane_only` via `--custom_mane_metadata_path` and `--custom_mane_pdb_dir`.
>
> - When invoking `oncodrive3d build-datasets`, supply:
>   - `--custom_mane_pdb_dir`: use this to provide your own predicted PDB structures (e.g., from [nf-core/proteinfold](https://nf-co.re/proteinfold/1.1.1/) or `final_bundle/pdbs/` produced by `update_samplesheet_and_structures.py`).
>   - `--custom_mane_metadata_path`: path to the corresponding `samplesheet.csv` (e.g., `final_bundle/samplesheet.csv`).
>
> See the documentation of [MANE Preprocessing Toolkit](tools/preprocessing/README.md) for the full workflow to expand coverage of the MANE associated structures.



## Running 3D clustering Analysis

For in depth information on how to obtain the required input data and for comprehensive information about the output, please refer to the [Input and Output Documentation](https://github.com/bbglab/oncodrive3d/tree/master/docs/run_input_output.md) of the 3D clustering analysis.  

### Input

- **Mutations file** (`required`): It can be either:
   - **<input_maf>**: A Mutation Annotation Format (MAF) file annotated with consequences (e.g., by using [Ensembl Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html)).
   - **<input_vep>**: The unfiltered output of VEP including annotations for all possible transcripts.

- **<mut_profile>** (`optional`): Dictionary including the normalized frequencies of mutations (*values*) in every possible trinucleotide context (*keys*), such as 'ACA>A', 'ACC>A', and so on.

---

> [!NOTE] 
> Examples of the input files are available in the [Test Input Folder](https://github.com/bbglab/oncodrive3d/tree/master/test/input).  
Please refer to these examples to understand the expected format and structure of the input files.

---

---

> [!NOTE]
> Oncodrive3D uses the mutational profile of the cohort to build an accurate background model. However, it’s not strictly required. If the mutational profile is not provided, the tool will use a simple uniform distribution as the background model for simulating mutations and scoring potential 3D clusters.

---

### Main Output

- **Gene-level output**: CSV file (`\<cohort>.3d_clustering_genes.csv`) containing the results of the analysis at the gene level. Each row represents a gene, sorted from the most significant to the least significant based on the 3D clustering analysis. The table also includes genes that were not analyzed, with the reason for exclusion provided in the `status` column.
  
- **Residue-level output**: CSV file (`<cohort>.3d_clustering_pos.csv`) containing the results of the analysis at the level of mutated residues. Each row corresponds to a mutated position within a gene and includes detailed information for each potential mutational cluster.


### Usage

```
Usage: oncodrive3d run [OPTIONS]

Examples:
  Basic run:
    oncodrive3d run -i <input_maf> -p <mut_profile> -d <build_folder> -C <cohort_name>
  
  Example of run using VEP output as input and MANE Select transcripts:
    oncodrive3d run -i <input_vep> -p <mut_profile> -d <build_folder> -C <cohort_name> \
                    --o3d_transcripts --use_input_symbols --mane

Options:
  -i, --input_path PATH            Path to the input file (MAF or VEP output) containing the 
                                   annotated mutations for the cohort. [required]
  -p, --mut_profile_path PATH      Path to the JSON file specifying the cohort's mutational 
                                   profile (192 key-value pairs).
  -m, --mutability_config_path PATH
                                   Path to the JSON file with information on mutability, used when
                                   mutation calls originate from cohorts with highly heterogeneous
                                   sequencing depth across sites.
  -o, --output_dir PATH            Path to the output directory for results. 
                                   Default: ./output/
  -d, --data_dir PATH              Path to the directory containing the datasets built in the 
                                   building datasets step. 
                                   Default: ./datasets/
  -n, --n_iterations INT           Number of densities to be simulated.
                                   Default: 10000
  -a, --alpha FLOAT                Significance threshold for residue and gene p-values.
                                   Default: 0.01
  -c, --cores INT                  Number of CPU cores to use. 
                                   Default: All available CPU cores
  -s, --seed INT                   Random seed for reproducibility.
  -v, --verbose                    Enables verbose output.
  -t, --cancer_type STR            Cancer type to include as metadata in the output file.
  -C, --cohort STR                 Cohort name for metadata and output file naming. 
  --no_fragments                   Disable processing of fragmented (AF-F) proteins.
  --only_processed                 Include only processed genes in the output.
  --thr_mapping_issue FLOAT        Threshold for the ratio of mutations with mapping issues 
                                   (out of structure, WT mismatch, zero mutability). A value 
                                   of 1 disables filtering mismatching mutations.
                                   Default: 0.1
  -P, --cmap_prob_thr FLOAT        Threshold for defining residues contacts based on distance 
                                   on predicted structure and predicted aligned error (PAE). 
                                   Default: 0.5
  --mane                           Prioritizes MANE Select transcripts when multiple 
                                   structures map to the same gene symbol.
  --o3d_transcripts                Filters mutations including only transcripts in Oncodrive3D 
                                   built datasets (requires VEP output as input file).
  --use_input_symbols              Update HUGO symbols in Oncodrive3D built datasets using the 
                                   input file's entries (requires VEP output as input file).
  -h, --help                       Show this message and exit.  
```

---

> [!NOTE]
> To maximize the number of matching transcripts between the input mutations and the AlphaFold predicted structures used by Oncodrive3D, it is recommended to use the unfiltered output of VEP (including all possible transcripts) as input, along with the flags `--o3d_transcripts` `--use_input_symbols` in the `oncodrive3d run` command.

### Handling Heterogeneous Sequencing Depth

Oncodrive3D can ingest **site-specific mutability tables** when a single mutational profile is not representative of the cohort (e.g., mutation calling performed on highly heterogeneous-depth datasets such as ultra depth **Duplex sequencing** panels commonly used in normal tissue analysis). Provide an indexed TSV describing per-site mutability together with a JSON config via `--mutability_config_path`. The run automatically switches from trinucleotide rates to per-position probabilities and tracks additional diagnostics (`Mut_zero_mut_prob`, `Pos_zero_mut_prob`, status `Mut_with_zero_prob/No_mutability`).  
> [!NOTE] 
> Please see the [Mutability-aware runs guide](docs/mutability.md) for the expected file formats, config schema, and troubleshooting tips.


### Building Annotations & Plotting

High-quality plots and enriched tables require a one-time `oncodrive3d build-annotations` step that retrieves AlphaFold structural features (secondary structure, surface accessibility), Pfam domains, UniProt features, and optional stability change predictions. Afterwards, `oncodrive3d plot` combines the run outputs with those annotations to create cohort summaries, per-gene multi-track panels, annotated CSVs, and association analyses between clusters and features.  
> [!NOTE] 
> Please refer to [Annotation & plotting workflow](docs/annotations_plotting.md) for prerequisites, command examples, and output descriptions.

---

### Running With Singularity

```
singularity pull oncodrive3d.sif docker://bbglab/oncodrive3d:latest
singularity exec oncodrive3d.sif oncodrive3d run -i <input_maf> -p <mut_profile> \ 
                                                 -d <build_folder> -C <cohort_name>
```


### Testing

To verify that Oncodrive3D is installed and configured correctly, you can perform a test run using the provided test input files: 

```
oncodrive3d run -d <build_folder> \
                -i ./test/input/maf/TCGA_WXS_ACC.in.maf \ 
                -p ./test/input/mut_profile/TCGA_WXS_ACC.sig.json \
                -o ./test/output/ -C TCGA_WXS_ACC
```

Check the output in the `test/output/` directory to ensure the analysis completes successfully.


## Parallel Processing on Multiple Cohorts

This repository provides a [Nextflow](https://www.nextflow.io/) pipeline to run Oncodrive3D in parallel across multiple cohorts, enabling efficient, reproducible and scalable analysis across datasets.  

For more information, refer to the [Oncodrive3D Pipeline](https://github.com/bbglab/oncodrive3d/tree/master/oncodrive3d_pipeline/) documentation.

### Usage

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
                                  `maf/` or `vep/` and `mut_profile/`. 
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
  --vep_input BOOL                Use `vep/` subdir as input and select transcripts matching 
                                  the Ensembl transcript IDs in Oncodrive3D built datasets. 
                                  Default: false
  --mane BOOL                     Prioritize structures corresponding to MANE transcrips if 
                                  multiple structures are associated to the same gene.
                                  Default: false
  --seed INT:                     Seed value for reproducibility.
                                  Default: 128
```


## License

Oncodrive3D is available to the general public subject to certain conditions described in its [LICENSE](LICENSE).


## Citation 
If you use Oncodrive3D in your research, please cite the original paper: [Oncodrive3D: fast and accurate detection of structural clusters of somatic mutations under positive selection](https://academic.oup.com/nar/article/53/15/gkaf776/8234003).
