// nextflow run main.nf --indir /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer/ --outdir /workspace/projects/clustering_3d/o3d_analysys/datasets/output/cancer/o3d_output/run_backtr_seq/ --cohort_pattern TCGA* --data_dir /workspace/projects/clustering_3d/clustering_3d/datasets_normal/
// nextflow run main.nf --indir /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer/ --cohort_pattern TCGA* --data_dir /workspace/nobackup/scratch/oncodrive3d/datasets
// nextflow run main.nf --indir /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer/ --outdir /workspace/projects/clustering_3d/o3d_analysys/datasets/output/cancer/o3d_output --data_dir /workspace/nobackup/scratch/oncodrive3d/datasets -profile conda
// nextflow run main.nf --indir /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer/ --outdir /workspace/projects/clustering_3d/o3d_analysys/datasets/output/cancer/o3d_output/no_mane --data_dir /workspace/nobackup/scratch/oncodrive3d/datasets_last_real -profile conda --cohort_pattern TCGA_WXS_BLCA

// Run no-MANE old input
// Run MANE new input -> add vep input 

// Add condition to 
// - use vep output as input
// - use o3d_transcripts 

input_files = params.vep_input ?
    "${params.indir}/{vep,mut_profile}/${params.cohort_pattern}{.vep.tsv.gz,.mutrate.json}" :
    "${params.indir}/{maf,mut_profile}/${params.cohort_pattern}{.in.maf,.mutrate.json}"
params.o3d_transcripts = params.vep_input ? true : params.o3d_transcripts
outdir = "${params.outdir}/${params.outsubdir}"  
 
log.info """\

    O n c o d r i v e - 3 D 
    =======================
    Input dir        : ${params.indir}
    Input files      : ${input_files}                      ## To del?
    Cohort pattern   : ${params.cohort_pattern}
    Outdir           : ${outdir}
    Datasets         : ${params.data_dir}
    Annotations      : ${params.annotations_dir}
    CPU cores        : ${params.cores}
    Memory           : ${params.memory}
    Max running      : ${params.max_running}
    Use VEP as input : ${params.vep_input}
    O3D transcripts  : ${params.o3d_transcripts}
    Generate plots   : ${params.plot}
    Seed             : ${params.seed}
    Container        : ${params.container}
    Profile          : ${workflow.profile}

    """
    .stripIndent()

process O3D_run {
    tag "O3D $cohort"
    // label 'process_high'
    debug true

    //errorStrategy 'ignore'        
    container params.container               
    cpus params.cores
    memory params.memory
    maxForks params.max_running
    publishDir outdir, mode:'copy'  //, pattern: "{$}"

    input:
    tuple val(cohort), path(inputs)

    output:
    tuple val(cohort), path("**genes.tsv"), path("**pos.tsv") , emit : o3d_result
    path("**.log")                                            , emit : log

    script:
    """
    oncodrive3D run -i ${inputs[0]} -p ${inputs[1]} -d ${params.data_dir} -C ${cohort} -o ${cohort} -s ${params.seed} -c ${params.cores}
    """
}

process O3D_plot {
    tag "Plot $cohort"
    //label 'process_low'
    debug true

    container params.container               
    cpus 4
    memory "10G"
    maxForks params.max_running

    input:
    tuple val(cohort), path(inputs), path(genes_tsv), path(pos_tsv)

    output:
    val(cohort)

    script:
    """
    oncodrive3D plot --output_tsv --non_significant -r $cohort -g $genes_tsv -p $pos_tsv -i ${inputs[0]} -m ${inputs[1]} -o $outdir/$cohort -d ${params.data_dir} -a ${params.annotations_dir}
    """
}

workflow {
    Channel
        .fromFilePairs(input_files, checkIfExists: true)
        .set { file_pairs_ch }

    O3D_run(file_pairs_ch)
        .set { run_ch }

    if (params.plot) { 
        file_pairs_ch
            .join(run_ch.o3d_result)
            .set { plot_ch }
        O3D_plot(plot_ch)
    }
}

workflow.onComplete {
    log.info ( workflow.success ? "\n3D-clustering analysis completed! --> $outdir/\n" : "FAILED!" )
}
