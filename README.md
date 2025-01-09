# Oncodrive3D 

**Oncodrive3D** is a computational method for analyzing patterns of somatic mutations across tumors. It identifies **three-dimensional (3D) clusters** of missense mutations and detects genes under **positive selection**.

---

## Installation

Install via PyPI:


```bash
pip install oncodrive3d
```

Install from source:

```bash
git clone https://github.com/bbglab/clustering_3d.git     # >>> Modify to oncodrive3D
cd clustering_3d                                          # >>> Modify to oncodrive3D
pip install .
```

## Building Datasets

This step is required after installation or whenever you need to generate datasets for a different organism or apply a specific threshold to define amino acid contacts.

Basic build:

```bash
oncodrive3d build-datasets -o <build_folder>
```

Add the `--mane` flag to include AlphaFold-predicted structures corresponding to MANE Select transcripts:

```bash
oncodrive3d build-datasets -o <build_folder> --mane
```

### Command Line Options:

- **-o, --output_dir <path>**: Path to the directory where the output files will be saved. *default:* `./datasets/`.

- **-s, --organism <str>**: Specifies the organism (`human` or `mouse`). *default:* `human`.

- **-m, --mane <flag: set to enable>**: Use structures predicted from MANE Select transcripts (applicable to Homo sapiens only).

- **-d, --distance_threshold <int>**: Distance threshold in Ångströms for defining amino acid contacts. *default:* `10`.

- **-c, --cores <int>**: Number of CPU cores for computation. *default:* All available CPU cores.

- **-v, --verbose <flag: set to enable>**: Enables verbose output.


## Running 3D-clustering Analysis

### Input

- **<input_mutations>** (`required`): Mutation Annotation Format (MAF) file annotated with consequences (e.g., by using [Ensembl Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html)).

- **<mut_profile>** (`optional`): Dictionary including the normalized frequencies of mutations (*values*) in every possible trinucleotide context (*keys*), such as 'ACA>A', 'ACC>A', and so on.

### Output

- **<cohort>.3d_clustering_genes.csv**: This is a Comma-Separated Values (CSV) file containing the results of the analysis at the gene level.
  
- **<cohort>.3d_clustering_pos.csv**: This is a Comma-Separated Values (CSV) file containing the results of the analysis at the level of mutated positions.


### Example Runs:

The `oncodrive3d run` command allows for flexible configurations based on the input data and desired options. Below are the typical scenarios:

Example of basic run:

```bash
oncodrive3d run -i <input_maf> -p <mut_profile> -d <build_folder> -C <cohort_name> -t <cancer_type>
```

Additional Options:

- Use MANE Select transcripts: Add the `--mane` flag to prioritize MANE transcripts.
- Use VEP output as input: Replace <input_maf> with <input_vep> and include the `--o3d_transcripts` and `--use_input_symbols` flags.
- Combine VEP input with MANE Select transcripts: Use both `--mane` and `--o3d_transcripts` flags together.

Example of run using VEP output as input and MANE Select transcripts:

```bash
oncodrive3d run -i <input_vep> -p <mut_profile> -d <build_folder> -C <cohort_name> -t <cancer_type> --o3d_transcripts --use_input_symbols --mane
```


#### Command Line Options:

- **-i, --input_path <path (required)>**: Path to the input file (MAF or VEP output) containing the annotated mutations for the cohort.

- **-p, --mut_profile_path <path>**: Path to the JSON file specifying the cohort's mutational profile (192 key-value pairs). If neither `--mut_profile_path` nor `--mutability_config_path` is provided, a uniform distribution across protein residues will be used to model neutral mutagenesis.

<!-- - **-m, --mutability_config_path <path>**: Path to the mutability configuration file, which contains integrated information about the mutation profile and sequencing depth of the cohor (required only for datasets generated with high variable depth between sites, such as those produced by Duplex Sequencing technology). -->

- **-o, --output_dir <path>**: Path to the output directory for results. *default:* `./output/`.

- **-d, --data_dir <path>**: Path to the directory containing the datasets built in the [building datasets](#building-datasets) step. *default:* `./datasets/`.

- **-c, --cores <int>**: Number of CPU cores to use. *default:* All available CPU cores.

- **-s, --seed <int>**: Random seed for reproducibility.

- **-v, --verbose <flag: set to enable>**: Enables verbose output.

- **-t, --cancer_type <str>**: Cancer type to include as metadata in the output file.

- **-C, --cohort <str>**: Cohort name for metadata and output file naming. 

- **-P, --cmap_prob_thr <float>**: Threshold for defining amino acid contacts based on distance on predicted structure and predicted aligned error (PAE). Default `0.5`.

- **--o3d_transcripts <flag: set to enable>**: Filters mutations to include only transcripts in Oncodrive3D built datasets (requires VEP output as input file).

- **--use_input_symbols <flag: set to enable>**: Update HUGO symbols in Oncodrive3D built datasets using the input file's entries (requires VEP output as input file).

- **--mane <flag: set to enable>**: Prioritizes MANE Select transcripts when multiple structures map to the same gene symbol.


### Running With Singularity

```bash
singularity pull oncodrive3d.sif docker://spellegrini87/oncodrive3d:1.0.2                                                              # TODO: TO UPDATE
singularity exec oncodrive3d.sif oncodrive3d run -i <input_maf> -p <mut_profile> -d <build_folder> -C <cohort_name>
```


## Testing

To verify that Oncodrive3D is installed and configured correctly, you can perform a test run using the provided test input files. 

Basic test run: 

```bash
   oncodrive3d run -d <build_folder> -i ./test/input/maf/TCGA_WXS_ACC.in.maf -p ./test/input/mut_profile/TCGA_WXS_ACC.sig.json -o ./test/output/ -C TCGA_WXS_ACC
```

Test run using VEP outout as Oncodrive3D input: 

```bash
   oncodrive3d run -d <build_folder> -i ./test/input/maf/TCGA_WXS_ACC.vep.tsv.gz -p ./test/input/mut_profile/TCGA_WXS_ACC.sig.json -o ./test/output/ -C TCGA_WXS_ACC --o3d_transcripts --use_input_symbols
```

Check the output in the `./test/output/` directory to ensure the analysis completes successfully.

The provided input files also serve as examples to demonstrate the expected structure and format of the input files. Users can refer to these files when preparing their own input data for analysis.


# Parallel Processing on Multiple Cohorts

Oncodrive3D can be run in parallel on multiple cohorts using [Nextflow](https://www.nextflow.io/). This approach enables efficient, reproducible and scalable analysis across datasets.

### Requirements

1. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) (version `23.04.3.5875` was used for testing).
2. Install either or both:
   - [Singularity](https://sylabs.io/guides/latest/user-guide/installation.html), or
   - [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).


### Setting Up and Running Oncodrive3D with Nextflow

#### Option 1: Using Singularity

Pull the Oncodrive3D Singularity images from Docker Hub.

```bash
   singularity pull build/containers/oncodrive3d.sif docker://spellegrini87/oncodrive3d:latest
```

#### Option 2: Using Conda

Ensure Oncodrive3D is installed in your Conda environment and update the `Profiles` section of the `nextflow.config` file to point to your Conda installation:

```groovy
    conda {
        process.executor = 'slurm'
        singularity.enabled = false
        conda.enabled = true
        process.conda = '<conda_environment_with_oncodrive3d>'
    }
```

Replace `<conda_environment_with_oncodrive3d>` with the path to your Conda environment.

#### Test Run

Run a test to ensure that everything is set up correctly and functioning as expected:

```bash
nextflow run main.nf -profile test,conda --data_dir <build_folder>
```

Replace `<build_folder>` with the path to the Oncodrive3D datasets built in the [building datasets](#building-datasets) step.
If you prefer to use Singularity, replace `conda` in the `-profile` argument with `container`.

#### Run on New Data

---

**⚠️ Note:**  
When using the Nextflow script, ensure that your input files are organized in the following directory structure:

```plaintext
input/
  ├── maf/
  │   └── <cohort>.in.maf
  ├── vep/
  │   └── <cohort>.vep.tsv.gz
  └── mut_profile/
      └── <cohort>.sig.json
```

- maf/: Contains mutation files with the `.in.maf` extension.
- vep/: Contains VEP annotation files with the `.vep.tsv.gz` extension, which include annotated mutations with all possible transcripts.
- mut_profile/: Contains mutational profile files with the `.sig.json` extension.

---

Example basic run using VEP output as input and MANE Select transcripts:

```bash
nextflow run main.nf -profile conda --data_dir <build_folder> --indir <input> --vep_input true --mane true
```

#### Command Line Options:

- **--indir <path>**: Path to the input directory including the subdirectories ``maf`` or ``vep`` and ``mut_profile``. *default:* ``${baseDir}/test/``

- **--outdir <path>**: Path to the output directory. *default:* ``run_<timestamp>/``

- **--cohort_pattern <str>**: Pattern expression to filter specific files within the input directory (e.g., 'TCGA*' would select only TCGA cohorts). *default:* ``*``

- **--data_dir <path>**: Path to the Oncodrive3D datasets directory, which includes the files compiled during the :ref:`building datasets` step. *default:* ``${baseDir}/datasets/``

- **--container <path>**: Path to the Singularity image with Oncodrive3D installation. *default:* ``${baseDir}/build/containers/oncodrive3d.sif``

- **--max_running <int>**: Maximum number of cohorts to process in parallel . *default:* ``5``

- **--cores <int>**: Number of CPU cores used to process each cohort. *default:* ``10``

- **--memory <str>**: Amount of memory allocated for processing each cohort. *default:* ``70GB``

- **--seed <int>**: Seed value for reproducibility. **default:** ``128``