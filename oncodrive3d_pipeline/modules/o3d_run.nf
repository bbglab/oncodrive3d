process O3D_RUN {
    tag "O3D $cohort"
    queue 'bigmem,normal'
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