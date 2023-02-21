#!/usr/bin/env nextflow

params.catalog = 'test_medicc_catalog.csv'

cohort_table = Channel.fromPath(params.catalog) |
    splitCsv(header:true)

input_file = cohort_table.map({row -> tuple(row.id, file(row.input_file, checkIfExists: true))})
medicc_args = cohort_table.map({row -> tuple(row.id, row.medicc_args)})
allele_specific = cohort_table.map({row -> tuple(row.id, row.allele_specific)})
output_directory = cohort_table.map({row -> tuple(row.id, row.output_directory)})

include { MEDICC_SIM } from './subworkflows/medicc_sim'

workflow {
    MEDICC_SIM(input_file, medicc_args, allele_specific, output_directory)
}
