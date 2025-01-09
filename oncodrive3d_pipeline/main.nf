#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { validatePaths } from './validation.nf'
include { ONCODRIVE3D } from './workflows/oncodrive3d'

workflow {
    validatePaths(params)
    ONCODRIVE3D()
}
