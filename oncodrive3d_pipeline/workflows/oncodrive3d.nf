include { O3D_RUN } from '../modules/o3d_run'
include { O3D_PLOT } from '../modules/o3d_plot'
include { O3D_CHIMERAX_PLOT } from '../modules/o3d_chimerax_plot'

workflow ONCODRIVE3D {
    main:
    // Define input channels
    input_files = params.vep_input ?
        "${params.indir}/{vep,mut_profile}/${params.cohort_pattern}{.vep.tsv.gz,.sig.json}" :
        "${params.indir}/{maf,mut_profile}/${params.cohort_pattern}{.in.maf,.sig.json}"

    log.info """
    O n c o d r i v e - 3 D 
    =======================
    Input dir          : ${params.indir}
    Input files        : ${input_files}                     
    Cohort pattern     : ${params.cohort_pattern}
    Output dir         : ${params.outdir}/${params.outsubdir}
    Datasets           : ${params.data_dir}
    Annotations        : ${params.annotations_dir}
    CPU cores          : ${params.cores}
    Memory             : ${params.memory}
    Max running        : ${params.max_running}
    Use VEP as input   : ${params.vep_input}
    MANE               : ${params.mane}
    Plots              : ${params.plot}
    ChimeraX plots     : ${params.chimerax_plot}
    Seed               : ${params.seed}
    Verbose            : ${params.verbose}
    Container          : ${params.container}
    Container ChimeraX : ${params.container_chimerax}
    Profile            : ${workflow.profile}
    """
    .stripIndent()
    
    Channel
        .fromFilePairs(input_files, checkIfExists: true)
        .map { cohort, files ->
            def mutFile = files.find { it.toString().endsWith(params.vep_input ? ".vep.tsv.gz" : ".in.maf") }
            def sigFile = files.find { it.toString().endsWith(".sig.json") }

            if (!mutFile || !sigFile) {
                error "Required files for cohort $cohort are missing: MUT file ($mutFile) or SIG file ($sigFile)."
            }

            return tuple(cohort, [mutFile, sigFile])
        }
        .set { file_pairs_ch }

    // Run O3D_RUN process
    O3D_RUN(file_pairs_ch)

    // Combine input files with O3D_RUN results for plotting
    plot_input_ch = file_pairs_ch.join(O3D_RUN.out.o3d_result)

    // Conditional execution of O3D_PLOT
    if (params.plot) {
        O3D_PLOT(plot_input_ch)
    }

    // Conditional execution of O3D_CHIMERAX_PLOT
    if (params.chimerax_plot) {
        O3D_CHIMERAX_PLOT(plot_input_ch)
    }

    emit:
        o3d_results = O3D_RUN.out.o3d_result
        o3d_logs = O3D_RUN.out.log

        plot_summary = params.plot ? O3D_PLOT.out.summary_plot : Channel.empty()
        plot_genes = params.plot ? O3D_PLOT.out.genes_plot : Channel.empty()
        plot_pos_annotated = params.plot ? O3D_PLOT.out.pos_annotated_csv : Channel.empty()
        plot_uniprot_feat = params.plot ? O3D_PLOT.out.uniprot_feat_csv : Channel.empty()
        plot_logs = params.plot ? O3D_PLOT.out.log : Channel.empty()

        plot_chimerax_defattr = params.chimerax_plot ? O3D_CHIMERAX_PLOT.out.chimerax_defattr : Channel.empty()
        plot_chimerax_plot = params.chimerax_plot ? O3D_CHIMERAX_PLOT.out.chimerax_plot : Channel.empty()
        plot_chimerax_log = params.chimerax_plot ? O3D_CHIMERAX_PLOT.out.log : Channel.empty()
}