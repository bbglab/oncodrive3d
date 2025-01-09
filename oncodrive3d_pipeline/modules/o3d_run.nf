process O3D_RUN {
    tag "O3D $cohort"
    
    container params.container
    cpus params.cores
    memory params.memory
    maxForks params.max_running
    publishDir "${params.outdir}/${params.outsubdir}", mode:'copy'

    // conda "bioconda::oncodrive3d"                                                                                    // TODO: Update and test
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?         // TODO: Update and test
    //     'https://depot.galaxyproject.org/singularity/oncodrive3d--py39hbf8eff0_0' :
    //     'quay.io/biocontainers/oncodrive3d--py39hbf8eff0_0' }"

    input:
    tuple val(cohort), path(inputs)

    output:
    tuple val(cohort), path("**genes.csv"), path("**pos.csv"), path("**mutations.processed.tsv"), path("**miss_prob.processed.json"), path("**seq_df.processed.tsv"), emit: o3d_result
    path "**.log", emit: log

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${cohort}"

    if (!params.data_dir) {
        error "Oncodrive3D built datasets path not provided. Please specify the path using --data_dir"
    }

    if (workflow.profile.contains('conda')) {
        if (params.data_dir == '/path/to/conda/environment/with/oncodrive3d') {
            error """
            Please update 'params.conda_env' in the nextflow.config file to the actual path of your Conda environment where Oncodrive3D is installed.
            Alternatively, you can provide the correct path as a command-line argument using:
            
            --conda_env <path_to_conda_environment>
            """
        }
    }

    """
    oncodrive3D run \\
        -i ${inputs[0]} \\
        -p ${inputs[1]} \\
        -d ${params.data_dir} \\
        -C ${cohort} \\
        -o ${prefix} \\
        -s ${params.seed} \\
        -c ${params.cores} \\
        ${params.ignore_mapping_issues ? '--thr_mapping_issue 1' : ''} \\
        ${params.verbose ? '-v' : ''} \\
        ${params.vep_input ? '--o3d_transcripts --use_input_symbols' : ''} \\
        ${params.mane ? '--mane' : ''} \\
        $args
    """
}