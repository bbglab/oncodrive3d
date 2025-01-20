# Oncodrive3D Nextflow Pipeline

This pipeline enables running Oncodrive3D in parallel across multiple cohorts using [Nextflow](https://www.nextflow.io/).

## Requirements

1. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) (version `23.04.3` was used for testing).
2. Install and set up either or both:
   - [Singularity](https://sylabs.io/guides/latest/user-guide/installation.html)  
      Pull the Oncodrive3D Singularity image from Docker Hub:

      ```
      singularity pull oncodrive3d.sif docker://bbglab/oncodrive3d:latest
      ```

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

Replace `<build_folder>` with the path to the Oncodrive3D datasets built in the [building datasets](#building-datasets) step.
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