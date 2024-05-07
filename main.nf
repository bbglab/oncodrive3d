// nextflow run main.nf --indir /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer/ --outdir /workspace/projects/clustering_3d/o3d_analysys/datasets/output/cancer/o3d_output/run_backtr_seq/ --cohort_pattern TCGA* --data_dir /workspace/projects/clustering_3d/clustering_3d/datasets_normal/
// nextflow run main.nf --indir /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer/ --cohort_pattern TCGA* --data_dir /workspace/nobackup/scratch/oncodrive3d/datasets
// nextflow run main.nf --indir /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer/ --outdir /workspace/projects/clustering_3d/o3d_analysys/datasets/output/cancer/o3d_output --data_dir /workspace/nobackup/scratch/oncodrive3d/datasets -profile conda
// nextflow run main.nf --indir /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer_202404/ --outdir /workspace/projects/clustering_3d/o3d_analysys/datasets/output/cancer/o3d_output/new_no_mane --data_dir /workspace/nobackup/scratch/oncodrive3d/datasets_last_real -profile conda --cohort_pattern PCAWG_WGS_ESO_ADENOCA --vep_input

// LAST HUMAN RAW nextflow run main.nf --indir /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer_202404/ --outdir /workspace/projects/clustering_3d/o3d_analysys/datasets/output/cancer_202404/o3d_output/human_raw --data_dir /workspace/nobackup/scratch/oncodrive3d/datasets_240506 -profile conda --vep_input true --verbose true
// LAST MANE RAW nextflow run main.nf --indir /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer_202404/ --outdir /workspace/projects/clustering_3d/o3d_analysys/datasets/output/cancer_202404/o3d_output/human_mane_raw --data_dir /workspace/nobackup/scratch/oncodrive3d/datasets_mane_240506 -profile conda --vep_input true --verbose true --mane true
// LAST HUMAN nextflow run main.nf --indir /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer_202404/ --outdir /workspace/projects/clustering_3d/o3d_analysys/datasets/output/cancer_202404/o3d_output/human --data_dir /workspace/nobackup/scratch/oncodrive3d/datasets_240506 -profile conda --verbose true
// LAST MANE nextflow run main.nf --indir /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer_202404/ --outdir /workspace/projects/clustering_3d/o3d_analysys/datasets/output/cancer_202404/o3d_output/human_mane --data_dir /workspace/nobackup/scratch/oncodrive3d/datasets_mane_240506 -profile conda --verbose true --mane true

// Example single run cancer tissue
// oncodrive3D run -i /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer_202404/vep/CBIOP_WXS_ANGS_UNTREAT_2020.vep.tsv.gz -p /workspace/projects/clustering_3d/o3d_analysys/datasets/input/cancer_202404/mut_profile/CBIOP_WXS_ANGS_UNTREAT_2020.sig.json -d /workspace/nobackup/scratch/oncodrive3d/datasets_240506 -C CBIOP_WXS_ANGS_UNTREAT_2020 -o /workspace/projects/clustering_3d/o3d_analysys/datasets/output/cancer_202404/o3d_output/human_raw/run_2024-05-07.CBIOP_WXS_ANGS_UNTREAT_2020 -s 26 -c 10 -v --o3d_transcripts --use_input_symbols


// Example single run normal tissue
// oncodrive3D run -i /workspace/datasets/transfer/ferriol_stefano/2024-04_data/single_sample/P19_0044_BDO_01.mutations.tsv -m /workspace/datasets/transfer/ferriol_stefano/2024-04_data/single_sample/oncodrive3d.mutability.conf -d /workspace/nobackup/scratch/oncodrive3d/datasets_240424 -C P19_0044_BDO_01 -o /workspace/projects/clustering_3d/o3d_analysys/datasets/output/normal/o3d_output/2024/P19_0044_BDO_01 -s 26 -c 20 -v
// oncodrive3D plot --output_tsv --non_significant -r $cohort -g $genes_tsv -p $pos_tsv -i ${inputs[0]} -m ${inputs[1]} -o $outdir/$cohort -d ${params.data_dir} -a ${params.annotations_dir} ${params.verbose ? '-v' : ''} 

// LAST PROCESSED oncodrive3D run -i /workspace/datasets/transfer/ferriol_stefano/2024-04_data/all_samples.mutations.tsv -m /workspace/datasets/transfer/ferriol_stefano/2024-04_data/oncodrive3d.mutability.conf -d /workspace/nobackup/scratch/oncodrive3d/datasets_240424 -C all_samples -o /workspace/projects/clustering_3d/o3d_analysys/datasets/output/normal/o3d_output/2024/all_samples_temp -s 26 -c 20 -v
// LAST RAW oncodrive3D run -i /workspace/datasets/transfer/ferriol_stefano/2024-04_data/all_samples.mutations.raw_vep.tsv -m /workspace/datasets/transfer/ferriol_stefano/2024-04_data/oncodrive3d.mutability.conf -d /workspace/nobackup/scratch/oncodrive3d/datasets_240424 -C all_samples_raw -o /workspace/projects/clustering_3d/o3d_analysys/datasets/output/normal/o3d_output/2024/all_samples_raw -s 26 -c 20 -v --o3d_transcripts --use_input_symbols

input_files = params.vep_input ?
    "${params.indir}/{vep,mut_profile}/${params.cohort_pattern}{.vep.tsv.gz,.sig.json}" :
    "${params.indir}/{maf,mut_profile}/${params.cohort_pattern}{.in.maf,.sig.json}"
outdir = "${params.outdir}/${params.outsubdir}"  
 
log.info """\

    O n c o d r i v e - 3 D 
    =======================
    Input dir        : ${params.indir}
    Input files      : ${input_files}                     
    Cohort pattern   : ${params.cohort_pattern}
    Outdir           : ${outdir}
    Datasets         : ${params.data_dir}
    Annotations      : ${params.annotations_dir}
    CPU cores        : ${params.cores}
    Memory           : ${params.memory}
    Max running      : ${params.max_running}
    Use VEP as input : ${params.vep_input}
    MANE             : ${params.mane}
    Generate plots   : ${params.plot}
    Seed             : ${params.seed}
    Verbose          : ${params.verbose}
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
    tuple val(cohort), path("**genes.csv"), path("**pos.csv"), path("**mutations.processed.tsv"), path("**seq_df.processed.tsv"), path("**miss_prob.processed.json"), emit : o3d_result
    path("**.log")                                                                                                                                                  , emit : log

    script:
    """
    oncodrive3D run -i ${inputs[0]} -p ${inputs[1]} -d ${params.data_dir} -C ${cohort} -o ${cohort} -s ${params.seed} -c ${params.cores} ${params.verbose ? '-v' : ''} ${params.vep_input ? '--o3d_transcripts --use_input_symbols' : ''} ${params.mane ? '--mane' : ''}
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
    oncodrive3D plot --output_tsv --non_significant -r $cohort -g $genes_tsv -p $pos_tsv -i ${inputs[0]} -m ${inputs[1]} -o $outdir/$cohort -d ${params.data_dir} -a ${params.annotations_dir} ${params.verbose ? '-v' : ''} 
    """
}


workflow {
    Channel
        .fromFilePairs(input_files, checkIfExists: true)
        .map { cohort, files ->
            // Find the VEP file as either .vep.tsv.gz or .in.maf based on params.vep_input
            def mutFile = files.find { it.toString().endsWith(params.vep_input ? ".vep.tsv.gz" : ".in.maf") }
            def sigFile = files.find { it.toString().endsWith(".sig.json") }

            // Ensure both files are present, throw an error if any is missing
            if (!mutFile || !sigFile) {
                throw new IllegalStateException("Required files for cohort $cohort are missing: MUT file ($mutFile) or SIG file ($sigFile).")
            }

            return tuple(cohort, [mutFile, sigFile]) // Return the explicitly ordered tuple
        }
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
