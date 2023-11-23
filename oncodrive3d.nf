params.indir = "${baseDir}/test"
params.cohort_pattern = "*"     
params.data_dir = "${baseDir}/datasets"            
params.container = "${baseDir}/build/containers/oncodrive3d.sif"
params.cores = 9
params.memory = "50G"
params.max_running = 5
params.seed = 128
input_files  = "${params.indir}/{maf,mut_profile}/${params.cohort_pattern}{.in.maf,.mutrate.json}"  
Date date = new Date()
String datePart = date.format("yyyyMMdd_mmss")
params.outdir = "run_${datePart}"
 
log.info """\

    O n c o d r i v e - 3 D 
    =======================
    Input dir     : ${params.indir}
    Cohort pattern: ${params.cohort_pattern}
    Outdir        : ${params.outdir}
    Datasets      : ${params.data_dir}
    CPU cores     : ${params.cores}
    Memory        : ${params.memory}
    Max running   : ${params.max_running}
    Seed          : ${params.seed}
    Container     : ${params.container}
    
    """
    .stripIndent()

process CLUSTERING {
    debug true

    //errorStrategy 'ignore'                       
    container params.container
    cpus params.cores
    memory params.memory
    maxForks params.max_running
    publishDir params.outdir, mode:'copy'
    tag "Clustering on $cohort"

    input:
    tuple val(cohort), path(inputs)

    output:
    path("${cohort}")
                                                                 /// REMEMBER TO CHANGE OPTION WITH NEW IMG        S -> s      u -> c
    script:
    """
    oncodrive3D run -i ${inputs[0]} -p ${inputs[1]} -d ${params.data_dir} -C ${cohort} -o ${cohort} -s ${params.seed} -c ${params.cores}
    """
}

workflow {
    Channel
        .fromFilePairs(input_files, checkIfExists: true)
        .set { file_pairs_ch }
    CLUSTERING(file_pairs_ch) 
}

workflow.onComplete {
    log.info ( workflow.success ? "\n3D-clustering analysis completed! --> $params.outdir/\n" : "FAILED!" )
}
