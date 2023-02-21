#!/usr/bin/env nextflow

params.catalog = 'test_medicc_ar_catalog.csv'

input_table = Channel.fromPath(params.catalog) |
    splitCsv(header:true)

tree = input_table.map({row -> tuple(row.id, file(row.tree_filename, checkIfExists: true))})
medicc_input = input_table.map({row -> tuple(row.id, file(row.medicc_input_file, checkIfExists: true))})
medicc_args = input_table.map({row -> tuple(row.id, row.medicc_args)})
allele_specific = input_table.map({row -> tuple(row.id, row.allele_specific)})
cell_list = input_table.map({row -> tuple(row.id, row.cell_list)})
output_directory = input_table.map({row -> tuple(row.id, row.output_directory)})

include { MEDICC_SITKA_SIM } from './subworkflows/medicc_sitka_sim'

workflow {
    MEDICC_SITKA_SIM(tree, medicc_input, medicc_args, allele_specific, cell_list, output_directory)
}
