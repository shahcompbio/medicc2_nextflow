#!/usr/bin/env nextflow

params.id ?: { params.id = 'output' }()

// Required arguments
params.tree_filename ?: { log.error "No tree provided. Make sure you have used the '--tree_filename' option."; exit 1 }()
tree = Channel.of(tuple(params.id, file(params.tree_filename, checkIfExists: true)))

params.signals_filename ?: { log.error "No signals data provided. Make sure you have used the '--signals_filename' option."; exit 1 }()
signals = Channel.of(tuple(params.id, file(params.signals_filename, checkIfExists: true)))

params.output_directory ?: { log.error "No copy number data provided. Make sure you have used the '--output_directory' option."; exit 1 }()
output_directory = Channel.of(tuple(params.id, params.output_directory))

params.segments_filename ?: { log.error "No segments provided. Make sure you have used the '--segments_filename' option."; exit 1 }()
segments = Channel.of(tuple(params.id, params.segments_filename))

// Optional arguments
params.cell_list = 'None'
cell_list = Channel.of(tuple(params.id, params.cell_list))

// Conditional arguments
params.allele_specific = ""
if (params.allele_specific == true){
    allele_specific = Channel.of(tuple(params.id, "--allele_specific")) // need to pass this flag to create_medicc_input
    medicc_args = Channel.of(tuple(params.id, ""))

} else {
    allele_specific = Channel.of(tuple(params.id, params.allele_specific))
    medicc_args = Channel.of(tuple(params.id, "--total-copy-numbers --input-allele-columns cn"))
}

include { MEDICC_SITKA } from './subworkflows/medicc_sitka'

workflow {
    MEDICC_SITKA(tree, signals, segments, medicc_args, allele_specific, cell_list, output_directory)
}
