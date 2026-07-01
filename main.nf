#!/usr/bin/env nextflow

/*
 * Root entry point required by Seqera Platform, which expects a repo-root
 * main.nf to recognise the repository as a Nextflow pipeline. The real
 * realbench pipeline lives under test/nextflow/ (kept there so the repo root
 * stays uncluttered); this thin wrapper imports and runs it so the pipeline
 * can be launched directly from the repo URL on Seqera / Tower.
 *
 * Locally you can still run the pipeline either way:
 *   nextflow run .                       -profile test   # via this wrapper
 *   nextflow run test/nextflow/main.nf   -profile test   # directly
 */

nextflow.enable.dsl = 2

include { REALBENCH } from './test/nextflow/main.nf'

workflow {
    REALBENCH()
}
