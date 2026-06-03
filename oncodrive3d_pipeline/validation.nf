// validation.nf

def validatePaths(params) {

    // Validate Conda environment
    if (workflow.profile.contains('conda')) {
        def condaPath = file(params.conda_env)
        if (!condaPath.exists() || !condaPath.isDirectory()) {
            error """
            \u001B[31mERROR: The specified Conda environment path does not exist or is not a directory:
            ${params.conda_env}

            Please ensure the correct path to the Conda environment containing Oncodrive3D is specified. 
            You can update 'params.conda_env' in the 'nextflow.config' file or provide it as a command-line argument:
            
            nextflow run main.nf --conda_env <path_to_conda_environment>\u001B[0m
            """
        }
    }

    // Validate input directory and its expected sub-directories
    if (!params.indir) {
        error """
        \u001B[31mERROR: Missing required parameter --indir.

        Provide the input directory containing the mutation files and mutational
        profiles, organised as:

            <indir>/${params.vep_input ? 'vep' : 'maf'}/<cohort>${params.vep_input ? '.vep.tsv.gz' : '.in.maf'}
            <indir>/mut_profile/<cohort>.sig.json

        nextflow run main.nf --indir <input_folder>\u001B[0m
        """
    }

    def inPath = file(params.indir)
    if (!inPath.exists() || !inPath.isDirectory()) {
        error """
        \u001B[31mERROR: The specified input directory does not exist or is not a directory:
        ${params.indir}\u001B[0m
        """
    }

    def mutSubdir = params.vep_input ? 'vep' : 'maf'
    [mutSubdir, 'mut_profile'].each { sub ->
        def subPath = file("${params.indir}/${sub}")
        if (!subPath.exists() || !subPath.isDirectory()) {
            error """
            \u001B[31mERROR: Expected sub-directory '${sub}/' not found under --indir:
            ${params.indir}/${sub}

            The input directory must contain '${mutSubdir}/' (mutation files) and
            'mut_profile/' (mutational profiles).\u001B[0m
            """
        }
    }

    // Validate Oncodrive3D datasets path
    def dataPath = file(params.data_dir)
    if (!dataPath.exists() || !dataPath.isDirectory()) {
        error """
        \u001B[31mERROR: The specified Oncodrive3D datasets path does not exist or is not a directory:
        ${params.data_dir}

        Please provide the path to the Oncodrive3D built datasets. 
        You can update 'params.data_dir' in the 'nextflow.config' file or provide it as a command-line argument: 

        nextflow run main.nf --data_dir <build_folder>\u001B[0m
        """
    }

    // Validate Oncodrive3D annotations datasets path
    if (params.plot == true) {
        def annotPath = file(params.annotations_dir)
        if (!annotPath.exists() || !annotPath.isDirectory()) {
            error """
            \u001B[31mERROR: The specified Oncodrive3D annotations datasets path does not exist or is not a directory:
            ${params.annotations_dir}

            Please provide the path to the Oncodrive3D built annotations datasets. 
            You can update 'params.annotations_dir' in the 'nextflow.config' file or provide it as a command-line argument: 

            nextflow run main.nf --annotations_dir <build_annotations_folder>\u001B[0m
            """
        }
    }
}