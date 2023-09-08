# Oncodrive3D 

__Oncodrive3D__ is a method designed to analyse patterns of somatic mutations 
across tumors to __identify three-dimensional (3D) clusters of missense mutations__ 
and __detect genes that are under positive selection__.

For detailed instructions on how to install, setup, and run the tool, how to obtain the required input data and for comprehensive information about the output, please refer to the [Oncodrive3D documentation](https://readthedocs.org).

## Installation

### Install using pip from PyPI

```bash
pip install oncodrive3d
```

### Install using pip from repository

```bash
git clone https://github.com/bbglab/clustering_3d.git
cd clustering_3d           # >>> Modify to oncodrive3D
pip install .
```

## Building datasets

```bash
oncodrive3D build-datasets
```

## Running 3D-clustering analysis

```bash
oncodrive3D run -i input.maf -p mut_profile.json -d build_folder/ -t cancer_type -C cohort_name
```

## Input & output

### Input

- **input.maf** (`required`): Mutation Annotation Format (MAF) file annotated with consequences (e.g., by using [Ensembl Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html)).

- **mut_profile.json** (`optional`): Dictionary including the normalized frequencies of mutations (`values`) in every possible trinucleotide context (`keys`), such as 'ACA>A', 'ACC>A', and so on.

### Output

- **cohort_filename.3d_clustering_genes.csv**: This is a Comma-Separated Values (CSV) file containing the results of the analysis at the gene level.
  
- **cohort_filename.3d_clustering_pos.csv**: This is a Comma-Separated Values (CSV) file containing the results of the analysis at the level of mutated positions.