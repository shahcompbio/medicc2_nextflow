#!/usr/bin/env nextflow

params.catalog = 'test_cohort_inputs.csv'

medicc_inputs = Channel.fromPath(params.catalog) |
    splitCsv(header:true)

tree = medicc_inputs.map({row -> tuple(row.id, file(row.tree_filename, checkIfExists: true))})
signals = medicc_inputs.map({row -> tuple(row.id, file(row.signals_filename, checkIfExists: true))})
segments = medicc_inputs.map({row -> tuple(row.id, file(row.segments_filename, checkIfExists: true))})
medicc_args = medicc_inputs.map({row -> tuple(row.id, row.medicc_args)})
allele_specific = medicc_inputs.map({row -> tuple(row.id, row.allele_specific)})
cell_list = medicc_inputs.map({row -> tuple(row.id, row.cell_list)})
output_directory = medicc_inputs.map({row -> tuple(row.id, row.output_directory)})

include { MEDICC_SITKA } from './subworkflows/medicc_sitka'

workflow {
    MEDICC_SITKA(tree, signals, segments, medicc_args, allele_specific, cell_list, output_directory)
}
