"""

#### EXAMPLE USAGE ####

python3 qmap_init.py -q submit.qmap -o /workspace/projects/clustering_3d/evaluation/tool_output/run_20230515_all_ext
python3 qmap_init.py -q submit.qmap -o /workspace/projects/clustering_3d/evaluation/tool_output/run_20230510 \
-e /home/odove/anaconda3/etc/profile.d/conda.sh

#######################

#### run qmap #########

qmap submit submit.qmap --max-running 30

#######################
"""


import argparse
import pandas as pd
import os


def init_submit_file(path_qmap_file, 
                     conda_sh_path, 
                     conda_env_name, 
                     memory = 10, cores = 1):
    """
    Create and initialize a submit.qmap file.
    """

    with open(path_qmap_file, "w") as file:
        file.write(f"[params]\nmemory = {memory}G\ncores = {cores}\n")
        file.write(f"\n[pre]\n")
        file.write(f'. "{conda_sh_path}"\n')
        file.write(f"conda activate {conda_env_name}\n\n[jobs]\n")


def add_job(script_dir, path_qmap_file, in_maf, in_mut_profile, output, cancer, cohort, num_iteration, ext_hits):
    """
    Add clustering_3d job to the qmap file.
    """

    command = f"python3 {script_dir}/main.py -i {in_maf} -o {output} -p {in_mut_profile} -H 0 -t {cancer} -C {cohort} -n {num_iteration} -e {ext_hits}"
    with open(path_qmap_file, "a") as file:
        file.write(command + "\n")
        

def init_parser():
    """
    Initialize parser for the main function.
    """

    SCRIPT_PATH = "/workspace/projects/clustering_3d/clustering_3d/scripts"
    COHORTS_PATH = "/workspace/projects/clustering_3d/evaluation/datasets/cohorts.tsv"
    IN_MAF = "/workspace/projects/clustering_3d/evaluation/datasets/input/maf"
    IN_MUT_PROF = "/workspace/projects/clustering_3d/evaluation/datasets/input/mut_profile"
    CONDA_SH = "/home/spellegrini/miniconda3/etc/profile.d/conda.sh"
    CONDA_ENV = "clustering_3d"

    parser = argparse.ArgumentParser()
    parser.add_argument("-q", "--qmap", help="Path to the qmap file to run jobs in parallel", type=str, required=True) 
    parser.add_argument("-o", "--output", help="Path to the output dir", type=str, required=True) 

    parser.add_argument("-s", "--script_dir", help="Path to dir including the scripts of the tool", type=str, default=SCRIPT_PATH)
    
    parser.add_argument("-e", "--conda_sh", help="Path to your own conda.sh file", type=str, default=CONDA_SH)
    parser.add_argument("-E", "--conda_env", help="Name of conda environment", type=str, default=CONDA_ENV) 

    parser.add_argument("-m", "--memory", help="GB of memory allocated to each job", type=int, default=10) 
    parser.add_argument("-c", "--cores", help="Number of cores allocated to each job", type=int, default=1) 

    parser.add_argument("-M", "--metadata", help="Path to the cohorts.tsv file inlcuding cohorts metadata", type=str, default=COHORTS_PATH) 
    parser.add_argument("-i", "--input_maf", help="Path to the input MAF file", type=str, default=IN_MAF) 
    parser.add_argument("-p", "--input_mut_profile", help="Path to the input mut_profile", type=str, default=IN_MUT_PROF) 
    
    parser.add_argument("-n", "--n_iterations", help = "Number of densities to be simulated", type=int, default=10000)
    parser.add_argument("-x", "--ext_hits",
                        help = "If 1 extend clusters to all mutated residues in the significant volumes, if 0 extend only to the ones having an anomaly > expected", 
                        type=int, default=1)

    return parser.parse_args()


def main():

    # Parser
    args = init_parser()
    qmap_file = args.qmap
    output = args.output
    script_dir = args.script_dir
    conda_sh = args.conda_sh
    conda_env = args.conda_env
    memory = args.memory
    cores = args.cores
    metadata = args.metadata
    input_maf = args.input_maf
    input_mut_profile = args.input_mut_profile
    num_iteration = args.n_iterations
    ext_hits = args.ext_hits

    # Create output folder if needed
    if not os.path.exists(output):
        os.makedirs(output)
    
    # Create a submit.qmap file
    init_submit_file(qmap_file, conda_sh, conda_env, memory, cores)

    # Add to qmap a 3D clustering job for each cohort with available MAF and mut_profile    
    metadata = pd.read_csv(metadata, sep="\t")
    i = 0
    for _, row in metadata.iterrows():

        tumor = row["CANCER_TYPE"]
        cohort = row["COHORT"]
        maf = f"{input_maf}/{cohort}.in.maf"
        mut_profile = f"{input_mut_profile}/{cohort}.mutrate.json"

        if os.path.isfile(maf) and os.path.isfile(mut_profile):
            add_job(script_dir, qmap_file, maf, mut_profile, output, tumor, cohort, num_iteration, ext_hits)
            i += 1

    print(f"{i} jobs added to {qmap_file}")


if __name__ == "__main__":
    main()


