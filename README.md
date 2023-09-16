# Oncodrive3D 

__Oncodrive3D__ is a method designed to analyse patterns of somatic mutations 
across tumors to identify __three-dimensional (3D) clusters__ of missense mutations 
and detect genes that are under __positive selection__.

For detailed instructions on how to install, setup, and run the tool, how to obtain the required input data and for comprehensive information about the output, please refer to the [Oncodrive3D documentation](https://readthedocs-toy.readthedocs.io/en/latest/).


## Installation

Install using pip from PyPI:

```bash
pip install oncodrive3d
```

Install using pip from repository:

```bash
git clone https://github.com/bbglab/clustering_3d.git
cd clustering_3d           # >>> Modify to oncodrive3D
pip install .
```

## Building datasets

This step is required once after installation or any time the user wants to 
build a different dataset (e.g., compile the datasets for a different organism 
or using a different threshold to define contacts between residues). 

```bash
oncodrive3D build-datasets -o build_folder/
```

- **-o, --output_dir <path>**: Specifies the build folder where files will be saved. Default: `datasets/`.

- **-s, --organism <str>**: Sets the organism name (`human` or `mouse`). Default: `human`.

- **-d, --distance_threshold <int>**: Distance threshold (Å) to define contact between amino acids. Default: `10`.

- **-c, --cores <int>**: Determines the number of CPU cores to use in the computation. Default: Number of available CPU cores.

- **-v, --verbose <flag: set to enable>**: Enables a more verbose output from the method.


## Running 3D-clustering analysis

```bash
oncodrive3D run -i input.maf -p mut_profile.json -d build_folder/ -t cancer_type -C cohort_name
```

- **-i, --input_maf_path <path (required)>**: Specifies the path to the MAF file of the cohort, including annotated mutations.

- **-p, --mut_profile_path <path>**: Specifies the path to the Mut profile of the cohort, which is a dictionary of 192 key-value pairs in JSON format.

- **-o, --output_dir <path>**: Sets the output directory. Default: `results/`.

- **-d, --data_dir <path>**: Sets the build folder, including the files compiled during the [building datasets](#building-datasets) step. Default: `datasets/`.

- **-c, --cores <int>**: Specifies the number of CPU cores to use in the computation. Default: Number of available CPU cores.

- **-s, --seed <int>**: Sets the seed to be used for reproducibility.

- **-v, --verbose <flag: set to enable>**: Enables a more verbose output from the method.

- **-t, --cancer_type <str>**: Specifies the cancer type used as metadata in the output file.

- **-C, --cohort <str>**: Specifies the cohort name used as metadata and filename for the output file.


## Input & output

### Input

- **input.maf** (`required`): Mutation Annotation Format (MAF) file annotated with consequences (e.g., by using [Ensembl Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html)).

- **mut_profile.json** (`optional`): Dictionary including the normalized frequencies of mutations (*values*) in every possible trinucleotide context (*keys*), such as 'ACA>A', 'ACC>A', and so on.

### Output

- **cohort_filename.3d_clustering_genes.csv**: This is a Comma-Separated Values (CSV) file containing the results of the analysis at the gene level.
  
- **cohort_filename.3d_clustering_pos.csv**: This is a Comma-Separated Values (CSV) file containing the results of the analysis at the level of mutated positions.


## Testing

To ensure that Oncodrive3D is correctly installed and configured, you can 
perform a test run using the provided test input files. 

```bash
   oncodrive3D run -i test/TCGA_WXS_ACC.in.maf -p test/TCGA_WXS_ACC.mutrate.json -o test/results/
```

## Parallel processing on multiple cohort

It is possible to run Oncodrive3D in parallel on multiple cohorts by using 
the provided [nextflow](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/#introduction) 
script:


```bash
   nextflow run oncodrive3d.nf --indir test/ --outdir test/results/
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

When using the nextflow script, it is important to ensure that your input 
*maf* and *mut profile* files are located in the same folder, as shown in 
``test/``. These files should have the extensions ``.in.maf`` 
and ``.mutrate.json``, respectively.