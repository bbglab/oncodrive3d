process O3D_CHIMERAX_PLOT {
    tag "ChimeraX plot $cohort"
    label 'process_low'
    queue 'bigmem,normal'
    maxForks params.max_running
    publishDir "${params.outdir}/${params.outsubdir}", mode:'copy'
    container 'docker.io/spellegrini87/oncodrive3d_chimerax:latest'

    input:
    tuple val(cohort), path(inputs), path(genes_csv), path(pos_csv), path(mutations_csv), path(miss_prob_json), path(seq_df_tsv)

    output:
    tuple val(cohort), path("**.chimerax/attributes/**.defattr"), optional: true, emit: chimerax_defattr
    tuple val(cohort), path("**.chimerax/plots/**.png"), optional: true, emit: chimerax_plot
    path "**.log", emit: log

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${cohort}"
    """
    oncodrive3D chimerax-plot \\
        -o $prefix \\
        -g $genes_csv \\
        -p $pos_csv \\
        -d ${params.data_dir} \\
        -s $seq_df_tsv \\
        -c $cohort \\
        --fragmented_proteins \\
        --transparent_bg \\
        ${params.verbose ? '-v' : ''} \\
        $args
    """
}