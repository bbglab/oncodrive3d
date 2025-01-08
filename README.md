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

## Building datasets

This step is required after installation or whenever you need to generate datasets for a different organism or apply a specific threshold to define amino acid contacts.

Basic build:

```bash
oncodrive3d build-datasets -o <build_folder>
```

Build with MANE Select transcripts:

```bash
oncodrive3d build-datasets -o <build_folder> --mane
```

### Command Line Options:

- **-o, --output_dir <path>**: Path to the directory where the output files will be saved. Default: `./datasets/`.

- **-s, --organism <str>**: Specifies the organism (`human` or `mouse`). Default: `human`.

- **-m, --mane <flag: set to enable>**: Use structures predicted from MANE Select transcripts (applicable to Homo sapiens only).

- **-d, --distance_threshold <int>**: Distance threshold in Ångströms for defining amino acid contacts. Default: `10`.

- **-c, --cores <int>**: Number of CPU cores for computation. Default: All available CPU cores.

- **-v, --verbose <flag: set to enable>**: Enables verbose output.


## Running 3D-clustering analysis

Basic run:

```bash
oncodrive3d run -i <input_maf> -p <mut_profile> -d <build_folder> -C <cohort_name>
```

Run with MANE Select transcripts:

```bash
oncodrive3d run -i <input_maf> -p <mut_profile> -d <build_folder> -C <cohort_name> --mane
```

Run using VEP outout as Oncodrive3D input:

```bash
oncodrive3d run -i <input_vep> -p <mut_profile> -d <build_folder> -C <cohort_name> --o3d_transcripts --use_input_symbols
```

Run using VEP outout as Oncodrive3D input and MANE Select transcripts:

```bash
oncodrive3d run -i <input_vep> -p <mut_profile> -d <build_folder> -C <cohort_name> --o3d_transcripts --use_input_symbols --mane
```

### Command Line Options:

- **-i, --input_path <path (required)>**: Path to the input file (MAF or VEP output) containing the annotated mutations for the cohort.

- **-p, --mut_profile_path <path>**: Path to the JSON file specifying the cohort's mutational profile (192 key-value pairs). If neither `--mut_profile_path` nor `--mutability_config_path` is provided, a uniform distribution across protein residues will be used to model neutral mutagenesis.

- **-m, --mutability_config_path <path>**Path to the mutability configuration file, which contains integrated information about the mutation profile and sequencing depth of the cohor (required only for datasets generated with high variable depth between sites, such as those produced by Duplex Sequencing technology).

- **-o, --output_dir <path>**: Output directory for results. Default: `./results/`.

- **-d, --data_dir <path>**: Path to the directory containing the datasets built in the [building datasets](#building-datasets) step. Default: `./datasets/`.

- **-c, --cores <int>**: Number of CPU cores to use. Default: All available CPU cores.

- **-s, --seed <int>**: Random seed for reproducibility.

- **-v, --verbose <flag: set to enable>**: Enables verbose output.

- **-t, --cancer_type <str>**: Cancer type to include as metadata in the output file.

- **-C, --cohort <str>**: Cohort name for metadata and output file naming. 

- **-P, --cmap_prob_thr <float>**: Threshold for defining amino acid contacts based on distance on predicted structure and predicted aligned error (PAE). Default `0.5`.

- **--o3d_transcripts <flag: set to enable>**: Filters mutations to include only transcripts in Oncodrive3D built datasets (requires VEP output as input file).

- **--use_input_symbols <flag: set to enable>**: Update HUGO symbols in Oncodrive3D built datasets using the input file's entries (requires VEP output as input file).

- **--mane <flag: set to enable>**: Prioritizes MANE Select transcripts when multiple structures map to the same gene symbol.


### Running with Singularity

```bash
singularity pull oncodrive3d.sif docker://spellegrini87/oncodrive3d:1.0.2                                                              # TODO: TO UPDATE
singularity exec oncodrive3d.sif oncodrive3d run -i <input_maf> -p <mut_profile> -d <build_folder> -C <cohort_name>
```


## Input & output

### Input

- **input.maf** (`required`): Mutation Annotation Format (MAF) file annotated with consequences (e.g., by using [Ensembl Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html)).

- **mut_profile.json** (`optional`): Dictionary including the normalized frequencies of mutations (*values*) in every possible trinucleotide context (*keys*), such as 'ACA>A', 'ACC>A', and so on.

- **mutability_config_path.json** (`optional`): # TODO: provide explaination: Dictionary of dictionary having...

### Output

- **cohort_filename.3d_clustering_genes.csv**: This is a Comma-Separated Values (CSV) file containing the results of the analysis at the gene level.
  
- **cohort_filename.3d_clustering_pos.csv**: This is a Comma-Separated Values (CSV) file containing the results of the analysis at the level of mutated positions.


## Testing

To verify that Oncodrive3D is installed and configured correctly, you can perform a test run using the provided test input files. 

Basic test run: 

```bash
   oncodrive3d run -d <build_folder> -i ./test/input/maf/TCGA_WXS_ACC.in.maf -p ./test/input/mut_profile/TCGA_WXS_ACC.sig.json -o ./test/output/ -C TCGA_WXS_ACC
```

Test run using VEP outout as Oncodrive3D input: 

```bash
   oncodrive3d run -d <build_folder> -i ./test/input/maf/TCGA_WXS_ACC.in.maf -p ./test/input/mut_profile/TCGA_WXS_ACC.sig.json -o ./test/output/ -C TCGA_WXS_ACC --o3d_transcripts --use_input_symbols
```

Check the output in the `./test/output/` directory to ensure the analysis completes successfully.

The provided input files also serve as examples to demonstrate the expected structure and format of the input files. Users can refer to these files when preparing their own input data for analysis.


# Parallel processing on multiple cohort

It is possible to run Oncodrive3D in parallel on multiple cohorts by using [nextflow](https://www.nextflow.io/).

1. Install [nextflow](https://www.nextflow.io/docs/latest/getstarted.html) and [Singularity](https://www.nextflow.io/docs/latest/getstarted.html) (versions `23.04.3.5875` and `3.5.3` were used respectively).

2. Pull Oncodrive3D image from Docker Hub:

```bash
   singularity pull build/containers/oncodrive3d.sif docker://spellegrini87/oncodrive3d:1.0.2
```

3. Run Oncodrive3D in parallel on multiple cohorts by using the provided nextflow script. For example:

```bash
   nextflow run oncodrive3d.nf --indir test/input --outdir test/output/
```

The nextflow script takes the following arguments:

--indir <path>   Input directory including the subdirectories ``maf`` and ``mut_profile``. Default: ``${baseDir}/test/``

--outdir <path>   Output directory. Default: ``run_<timestamp>/``

--cohort_pattern <str>   Pattern expression to select specific files within the input directory (e.g., 'TCGA*' would select only TCGA cohorts). Default: ``*``

--data_dir <path>   Build folder including the files compiled during the :ref:`building datasets` step. Default: ``${baseDir}/datasets/``

--container <path>   Singularity image with installation of Oncodrive3D. Default: ``${baseDir}/build/containers/oncodrive3d.sif``

--max_running <int>   Maximum number of cohorts allowed to be processed in parallel . Default: ``5``

--cores <int>   CPU cores used to process each cohort. Default: ``9``

--memory <str>   Memory used to process each cohort. Default: ``50GB``

--seed <int>   Seed to be used for reproducibility. Default: ``128``

When using the nextflow script, it is important to ensure that your input *maf* and *mut profile* files are located in the same folder, as shown in 
``test/input``. These files should have the extensions ``.in.maf`` and ``.sig.json``, respectively.


## Quick interpretation of the analysis

You can generate plots for a quick interpretation of the 3D clustering analysis 
performed by Oncodrive-3D. The plots can be simple or annotated with structural 
and genomics features. To generate annotated plots, it is required (once) to 
build the annotations datasets.

### Installation of external software for annotations

Install PDB_Tool to determine solvent accessibility and secondary structures 
from the PDB files.

1. Clone PDB_Tool:

```bash
git clone https://github.com/realbigws/PDB_Tool
cd PDB_Tool/source_code
make
```

2. Open your configuration file for the Bash shell environment:

```bash
nano ~/.bashrc
```

3. Export the path for PDB_Tool to your enviroment variable by adding the 
following line to your configuraion file (change `/path/to/PDB_Tool` to the 
your actual path to PDB_Tool):

```bash
export PATH="$PATH:/path/to/PDB_Tool"
```

### Building annotations

This step is required once, only to enable Oncodrive3D to produce annotated plots.
It is not required to produce simple plot nor to run the 3D-clustering analysis.

```bash
oncodrive3D build-annotations -o annotation_folder/
```

- **-o, --output_dir <path>**: Specifies the annotation folder where files will be saved. Default: `./annotations/`.

- **-c, --cores <int>**: Determines the number of CPU cores to use in the computation. Default: Number of available CPU cores.

- **-v, --verbose <flag: set to enable>**: Enables a more verbose output from the method.


### Generating plots

```bash
oncodrive3D plot -i oncodrive3d_result/ -f filename
```
