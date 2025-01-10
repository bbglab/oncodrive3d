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

    // Validate Oncodrive3D Singularity image
    if (workflow.containerEngine == 'singularity') {
        def imagePath = file(params.container)
        if (!imagePath.exists()) {
            error """
            \u001B[31mERROR: The specified Oncodrive3D Singularity image does not exist:
            ${params.container}

            Please provide the path to the Oncodrive3D Singularity image. 
            You can update 'params.container' in the 'nextflow.config' file or provide it as a command-line argument: 

            nextflow run main.nf --container <singularity_image>\u001B[0m
            """
        }
    }

    // Validate Oncodrive3D with ChimeraX Singularity image
    if (params.chimerax_plot == true) {
        def image_chimeraxPath = file(params.container_chimerax)
        if (!image_chimeraxPath.exists() || !image_chimeraxPath.isDirectory()) {
            error """
            \u001B[31mERROR: The specified ChimeraX Singularity image does not exist:
            ${params.container_chimerax}

            Please provide the path to the ChimeraX Singularity image. 
            You can update 'params.container_chimerax' in the 'nextflow.config' file or provide it as a command-line argument: 

            nextflow run main.nf --container_chimerax <singularity_image>\u001B[0m
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