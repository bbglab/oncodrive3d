# Oncodrive3D 

**Oncodrive3D** is a fast and accurate computational method designed to analyze patterns of somatic mutation across tumors, with the goal of identifying **three-dimensional (3D) clusters** of missense mutations and detecting genes under **positive selection**. 

The method leverages **AlphaFold 2-predicted protein structures** and Predicted Aligned Error (PAE) to define residue contacts within the protein's 3D space. When available, it integrates **mutational profiles** to build an accurate background model of neutral mutagenesis. By applying a novel **rank-based statistical approach**, Oncodrive3D scores potential 3D clusters and computes empirical p-values."

---

## License

Oncodrive3D is available to the general public subject to certain conditions described in its [LICENSE](LICENSE).

## Installation

Install via PyPI:


```bash
pip install oncodrive3d
```

Alternatively, you can obtain the latest code from the repository and install it for development with pip:

```bash
git clone https://github.com/bbglab/oncodrive3d.git
cd oncodrive3d
python -m venv .venv
.venv/bin/pip install -e .
source .venv/bin/activate
oncodrive3d --help
```

> [!NOTE]
> If you install a modern build tool like [uv](https://github.com/astral-sh/uv),
> you can simply do this:
> ```bash
> git clone https://github.com/bbglab/clustering_3d.git
> cd clustering_3d
> uv run oncodrive3d --help
> ```

## Building Datasets

This step is required after installation or whenever you need to generate datasets for a different organism or apply a specific threshold to define amino acid contacts.

> [!NOTE]
> The first time that you run Oncodrive3D building dataset step with a given reference genome, it
> will download it from our servers. By default the downloaded datasets go to
> `~/.bgdata`. If you want to move these datasets to another folder you have to
> define the system environment variable `BGDATA_LOCAL` with an export command.

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
  -d, --distance_threshold INT    Distance threshold (Å) for defining residues contacts. 
                                  Default: 10
  -c, --cores INT                 Number of CPU cores for computation. 
                                  Default: All available CPU cores
  -v, --verbose                   Enables verbose output.
  -h, --help                      Show this message and exit.  
```

## Running 3D clustering Analysis

For in depth information on how to obtain the required input data and for comprehensive information about the output, please refer to the [Input and Output Documentation](https://github.com/bbglab/oncodrive3d/tree/master/docs/input_output.md).  

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

- **\<cohort>.3d_clustering_genes.csv**: A Comma-Separated Values (CSV) file containing the results of the analysis at the gene level. Each row represents a gene, sorted from the most significant to the least significant based on the 3D clustering analysis. The table also includes genes that were not analyzed, with the reason for exclusion provided in the `status` column.
  
- **\<cohort>.3d_clustering_pos.csv**: A Comma-Separated Values (CSV) file containing the results of the analysis at the level of mutated residues. Each row corresponds to a mutated position within a gene and includes detailed information for each mutational cluster.


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
  -o, --output_dir PATH            Path to the output directory for results. 
                                   Default: ./output/
  -d, --data_dir PATH              Path to the directory containing the datasets built in the 
                                   building datasets step. 
                                   Default: ./datasets/
  -c, --cores INT                  Number of CPU cores to use. 
                                   Default: All available CPU cores
  -s, --seed INT                   Random seed for reproducibility.
  -v, --verbose                    Enables verbose output.
  -t, --cancer_type STR            Cancer type to include as metadata in the output file.
  -C, --cohort STR                 Cohort name for metadata and output file naming. 
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

Oncodrive3D can be run in parallel on multiple cohorts using [Nextflow](https://www.nextflow.io/). This approach enables efficient, reproducible and scalable analysis across datasets.

### Requirements

1. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) (version `23.04.3` was used for testing).
2. Install either or both:
   - [Singularity](https://sylabs.io/guides/latest/user-guide/installation.html)
   - [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)


#### Option 1: Using Singularity

Pull the Oncodrive3D Singularity images from Docker Hub.

```
singularity pull oncodrive3d.sif docker://bbglab/oncodrive3d:latest
```

#### Option 2: Using Conda

Ensure Oncodrive3D is installed in your Conda environment and update the `params` section of the `nextflow.config` file to point to your Conda installation:

```groovy
params {
    ...
    conda_env = '/path/to/conda/environment/with/oncodrive3d' 
    ...
}
```

Replace `/path/to/conda/environment/with/oncodrive3d` with the path to your Conda environment. Alternatively, you can provide it as a command-line argument.

### Test Run

Run a test to ensure that everything is set up correctly and functioning as expected:

```
cd oncodrive3d_pipeline
nextflow run main.nf -profile test,container --data_dir <build_folder>
```

Replace `<build_folder>` with the path to the Oncodrive3D datasets built in the [building datasets](#building-datasets) step.
If you prefer to use Conda, replace `container` in the `-profile` argument with `conda`.

### Usage

---

> [!WARNING]
> When using the Nextflow script, ensure that your input files are organized in the following directory structure:
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
  Basic run:
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
  --container PATH                Path to the Singularity image with Oncodrive3D installation. 
                                  Default: ${baseDir}/../oncodrive3d.sif
  --max_running INT               Maximum number of cohorts to process in parallel.
                                  Default: 5
  --cores INT                     Number of CPU cores used to process each cohort. 
                                  Default: 10
  --memory STR                    Amount of memory allocated for processing each cohort. 
                                  Default: 70GB
  --seed INT:                     Seed value for reproducibility.
                                  Default: 128
```