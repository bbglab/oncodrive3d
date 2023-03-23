"""
#### run qmap ####


python3 qmap_init.py -q submit.qmap -o /workspace/projects/clustering_3d/evaluation/tool_output/run_23.03.23
qmap submit submit.qmap --max-running 30

##################
"""


import argparse
import pandas as pd
import os


def add_job(path_qmap_file, in_maf, in_mut_profile, output, cancer, cohort):
    """
    Add clustering_3d job to the qmap file.
    """

    command = f"python3 main.py -i {in_maf} -o {output} -p {in_mut_profile} -H 0 -t {cancer} -C {cohort}"
    with open(path_qmap_file, "a") as file:
        file.write(command + "\n")


def main():

    COHORTS_PATH = "/workspace/projects/clustering_3d/evaluation/datasets/cohorts.tsv"
    IN_MAF = "/workspace/projects/clustering_3d/evaluation/datasets/input/maf"
    IN_MUT_PROF = "/workspace/projects/clustering_3d/evaluation/datasets/input/mut_profile"

    # Parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-q", "--qmap", help="Path to the qmap file to run jobs in parallel", type=str, required=True) 
    parser.add_argument("-o", "--output", help="Path to the output dir", type=str, required=True) 
    parser.add_argument("-m", "--metadata", help="Path to the cohorts.tsv file inlcuding cohorts metadata", type=str, default=COHORTS_PATH) 
    parser.add_argument("-i", "--input_maf", help="Path to the input MAF file", type=str, default=IN_MAF) 
    parser.add_argument("-p", "--input_mut_profile", help="Path to the input mut_profile", type=str, default=IN_MUT_PROF) 

    args = parser.parse_args()
    qmap = args.qmap
    output = args.output
    metadata = args.metadata
    input_maf = args.input_maf
    input_mut_profile = args.input_mut_profile
    
    # Add to qmap a 3D clustering job for each cohort with available MAF and mut_profile
    metadata = pd.read_csv(COHORTS_PATH, sep="\t")
    for i, row in metadata.iterrows():

        tumor = row["CANCER_TYPE"]
        cohort = row["COHORT"]
        maf = f"{input_maf}/{cohort}.in.maf"
        mut_profile = f"{input_mut_profile}/{cohort}.mutrate.json"

        if os.path.isfile(maf) and os.path.isfile(mut_profile):
            add_job(qmap, maf, mut_profile, output, tumor, cohort)

    print(f"{i+1} jobs added to {qmap}")


if __name__ == "__main__":
    main()


