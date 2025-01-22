process O3D_PLOT {
    tag "Plot $cohort"
    label 'process_low'
    queue 'bigmem,normal'
    container params.container
    maxForks params.max_running
    publishDir "${params.outdir}/${params.outsubdir}", mode:'copy'

    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?       // TODO: Update and test
    //     'https://depot.galaxyproject.org/singularity/oncodrive3d--py39hbf8eff0_0' :
    //     'quay.io/biocontainers/oncodrive3d--py39hbf8eff0_0' }"

    input:
    tuple val(cohort), path(inputs), path(genes_csv), path(pos_csv), path(mutations_csv), path(miss_prob_json), path(seq_df_tsv)

    output:
    tuple val(cohort), path("**.summary_plot.png"), optional: true, emit: summary_plot
    tuple val(cohort), path("**.genes_plots/**.png"), optional: true, emit: genes_plot
    tuple val(cohort), path("**.associations_plots/**.logodds_plot.png"), optional: true, emit: logodds_plot
    tuple val(cohort), path("**.associations_plots/**.volcano_plot.png"), optional: true, emit: volcano_plot
    tuple val(cohort), path("**.associations_plots/**.volcano_plot_gene.png"), optional: true, emit: volcano_plot_gene
    tuple val(cohort), path("**.3d_clustering_pos.annotated.csv"), optional: true, emit: pos_annotated_csv
    tuple val(cohort), path("**.uniprot_feat.tsv"), optional: true, emit: uniprot_feat_tsv
    tuple val(cohort), path("**.logreg_result.tsv"), optional: true, emit: logreg_tsv
    path "**.log", emit: log

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${cohort}"
    """
    oncodrive3D plot \\
        -g $genes_csv \\
        -p $pos_csv \\
        -i $mutations_csv \\
        -m $miss_prob_json \\
        -s $seq_df_tsv \\
        -d ${params.data_dir} \\
        -a ${params.annotations_dir} \\
        -o $prefix \\
        -c $cohort \\
        --output_csv \\
        ${params.verbose ? '-v' : ''} \\
        $args
    """
}