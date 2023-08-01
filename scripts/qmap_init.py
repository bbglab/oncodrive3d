"""
qmap initializer

#### EXAMPLE USAGE ####

python3 qmap_init.py -q submit_pcmap_0.5.qmap -o /workspace/projects/clustering_3d/evaluation/tool_output/run_20230727_pcmap_0.5_final \
    -p /workspace/projects/clustering_3d/evaluation/datasets/input/mut_profile -c  \
        /workspace/projects/clustering_3d/clustering_3d/datasets_frag/prob_cmaps/  \
            -d /workspace/projects/clustering_3d/clustering_3d/datasets_frag/confidence.csv \
                -s /workspace/projects/clustering_3d/clustering_3d/datasets_frag/seq_for_mut_prob.csv -u 10 -m 55 -r 128 -P 0.5

#######################


#### run qmap #########

qmap submit submit_pcmap_0.5.qmap --max-running 10

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


def add_job(script_dir, path_qmap_file, 
            in_maf, in_mut_profile, output,
            seq_df, cmap_path, plddt_path, cores,
            cancer, cohort, num_iteration, seed,
            cmap_prob_thr, pae_path):
    """
    Add clustering_3d job to the qmap file.
    """

    if in_mut_profile is not None:
        flag_mut_profile = f"-p {in_mut_profile}"
    else:
        flag_mut_profile = ""
    if seed is not None:
        flag_seed = f"-S {seed}"
    else:
        flag_seed = ""
    command = f"python3 {script_dir}/main.py -i {in_maf} -o {output} {flag_mut_profile} -s {seq_df} -c {cmap_path} -d {plddt_path} -t {cancer} -C {cohort} -u {cores} -n {num_iteration} -P {cmap_prob_thr} {flag_seed} -e {pae_path}"
    with open(path_qmap_file, "a") as file:
        file.write(command + "\n")
        

def init_parser():
    """
    Initialize parser for the main function.
    """

    SCRIPT_PATH = "/workspace/projects/clustering_3d/clustering_3d/scripts/"
    COHORTS_PATH = "/workspace/projects/clustering_3d/evaluation/datasets/cohorts.tsv"
    IN_MAF = "/workspace/projects/clustering_3d/evaluation/datasets/input/maf"
    IN_SEQ = "/workspace/projects/clustering_3d/clustering_3d/datasets/seq_for_mut_prob.csv"
    IN_CMAP = "/workspace/projects/clustering_3d/clustering_3d/datasets/cmaps/"
    IN_PLDDT = "/workspace/projects/clustering_3d/clustering_3d/datasets/confidence.csv"
    
    CONDA_SH = "/home/spellegrini/miniconda3/etc/profile.d/conda.sh"
    CONDA_ENV = "clustering_3d"

    parser = argparse.ArgumentParser()
    parser.add_argument("-q", "--qmap", help="Path to the qmap file to run jobs in parallel", type=str, required=True) 
    parser.add_argument("-o", "--output", help="Path to the output dir", type=str, required=True) 
    
    parser.add_argument("-M", "--metadata", help="Path to the cohorts.tsv file inlcuding cohorts metadata", type=str, default=COHORTS_PATH) 
    parser.add_argument("-i", "--input_maf", help="Path to the input MAF file", type=str, default=IN_MAF) 
    parser.add_argument("-p", "--input_mut_profile", help="Path to the input mut_profile", type=str) 
    
    parser.add_argument("-s", "--seq_df", help = "Path to the dataframe including DNA and protein seq of all gene/proteins (all AF predicted ones)", type=str, default = IN_SEQ)       
    parser.add_argument("-c", "--cmap_path", help = "Path to the directory containting the contact map of each protein", type=str, default = IN_CMAP)
    parser.add_argument("-d", "--plddt_path", help = "Path to the pandas dataframe including the AF model confidence of all proteins", type=str, default = IN_PLDDT)
    parser.add_argument("-a", "--pae_path", help = "Path to the directory including the AF Predicted Aligned Error (PAE) .npy files", type=str)

    parser.add_argument("-S", "--script_dir", help="Path to dir including the scripts of the tool", type=str, default=SCRIPT_PATH)
    
    parser.add_argument("-e", "--conda_sh", help="Path to your own conda.sh file", type=str, default=CONDA_SH)
    parser.add_argument("-E", "--conda_env", help="Name of conda environment", type=str, default=CONDA_ENV) 

    parser.add_argument("-r", "--seed", help="Set seed to ensure reproducible results", type=int)
    parser.add_argument("-m", "--memory", help="GB of memory allocated to each job", type=int, default=10) 
    parser.add_argument("-u", "--cores", help="Number of cores allocated to each job", type=int, default=1) 
    
    parser.add_argument("-n", "--n_iterations", help = "Number of densities to be simulated", type=int, default=10000)
    parser.add_argument("-P", "--cmap_prob_thr", help = "Threshold to define AAs contacts based on distance on predicted structure and PAE", type=float, default=0.5)


    return parser.parse_args()


def main():

    # Parser
    args = init_parser()
    qmap_file = args.qmap
    output = args.output
    
    seq_df = args.seq_df
    cmap_path = args.cmap_path
    plddt_path = args.plddt_path
    
    script_dir = args.script_dir
    conda_sh = args.conda_sh
    conda_env = args.conda_env
    seed = args.seed
    memory = args.memory
    cores = args.cores
    metadata = args.metadata
    input_maf = args.input_maf
    input_mut_profile = args.input_mut_profile
    num_iteration = args.n_iterations
    cmap_prob_thr = args.cmap_prob_thr
    pae_path = args.pae_path

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
        
        if cohort.startswith("ICGC"):
            continue

        if input_mut_profile is not None:
            mut_profile = f"{input_mut_profile}/{cohort}.mutrate.json"
            if os.path.isfile(maf) and os.path.isfile(mut_profile):
                add_job(script_dir, qmap_file, 
                        maf, mut_profile, output,
                        seq_df, cmap_path, plddt_path, cores,
                        tumor, cohort, num_iteration, seed, 
                        cmap_prob_thr, pae_path)
                i += 1
        else:
            mut_profile = None
            if os.path.isfile(maf):
                add_job(script_dir, qmap_file, 
                        maf, mut_profile, output,
                        seq_df, cmap_path, plddt_path, cores,
                        tumor, cohort, num_iteration, seed, 
                        cmap_prob_thr, pae_path)
                i += 1
            
    print(f"{i} jobs added to {qmap_file}")


if __name__ == "__main__":
    main()


