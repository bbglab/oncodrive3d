/// nextflow run main.nf --indir /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer/ --outdir /workspace/projects/clustering_3d/o3d_analysys/datasets/output/cancer/o3d_output/run_backtr_seq/ --cohort_pattern TCGA* --data_dir /workspace/projects/clustering_3d/clustering_3d/datasets_normal/
/// nextflow run main.nf --indir /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer/ --cohort_pattern TCGA* --data_dir /workspace/nobackup/scratch/oncodrive3d/datasets
/// nextflow run main.nf --indir /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer/ --data_dir /workspace/nobackup/scratch/oncodrive3d/datasets --outdir /workspace/projects/clustering_3d/o3d_analysys/datasets/output/cancer/o3d_output/ -profile conda

input_files = "${params.indir}/{maf,mut_profile}/${params.cohort_pattern}{.in.maf,.mutrate.json}"
outdir = "${params.outdir}/${params.outsubdir}"  
 
log.info """\

    O n c o d r i v e - 3 D 
    =======================
    Input dir      : ${params.indir}
    Cohort pattern : ${params.cohort_pattern}
    Outdir         : ${outdir}
    Datasets       : ${params.data_dir}
    CPU cores      : ${params.cores}
    Memory         : ${params.memory}
    Max running    : ${params.max_running}
    Seed           : ${params.seed}
    Container      : ${params.container}
    Profile        : ${workflow.profile}

    """
    .stripIndent()

process O3D {
    debug true

    //errorStrategy 'ignore'        
    container params.container               
    cpus params.cores
    memory params.memory
    maxForks params.max_running
    publishDir outdir, mode:'copy'
    tag "Clustering on $cohort"

    input:
    tuple val(cohort), path(inputs)

    output:
    path("${cohort}")

    script:
    """
    oncodrive3D run -i ${inputs[0]} -p ${inputs[1]} -d ${params.data_dir} -C ${cohort} -o ${cohort} -s ${params.seed} -c ${params.cores}
    """
}

workflow {
    Channel
        .fromFilePairs(input_files, checkIfExists: true)
        .set { file_pairs_ch }
    O3D(file_pairs_ch) 
}

workflow.onComplete {
    log.info ( workflow.success ? "\n3D-clustering analysis completed! --> $outdir/\n" : "FAILED!" )
}
