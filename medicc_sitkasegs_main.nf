#!/usr/bin/env nextflow

params.id ?: { params.id = 'output' }()

// Required arguments
params.signals_filename ?: { log.error "No signals data provided. Make sure you have used the '--signals_filename' option."; exit 1 }()
signals = Channel.of(tuple(params.id, file(params.signals_filename, checkIfExists: true)))

params.segments_filename ?: { log.error "No segments provided. Make sure you have used the '--segments_filename' option."; exit 1 }()
segments = Channel.of(tuple(params.id, file(params.segments_filename, checkIfExists: true)))

params.output_directory ?: { log.error "No copy number data provided. Make sure you have used the '--output_directory' option."; exit 1 }()
output_directory = Channel.of(tuple(params.id, params.output_directory))

// Optional argument
params.medicc_args = """-j 24 --input-type t --verbose --plot none --no-plot-tree --events \
--chromosomes-bed /data1/shahs3/isabl_data_lake/software/dependencies/medicc2/medicc/objects/hg19_chromosome_arms.bed \
--regions-bed /data1/shahs3/isabl_data_lake/software/dependencies/medicc2_nextflow/Davoli_2013_TSG_OG_genes_hg37.bed"""

params.cell_list = "None"
cell_list = Channel.of(tuple(params.id, params.cell_list))

// Conditional arguments
params.allele_specific = ""
if (params.allele_specific == true){
    allele_specific = Channel.of(tuple(params.id, "--allele_specific")) // need to pass this flag to create_medicc_input
    medicc_args = Channel.of(tuple(params.id, params.medicc_args))

} else {
    allele_specific = Channel.of(tuple(params.id, params.allele_specific))
    medicc_args = Channel.of(tuple(params.id, params.medicc_args + " --total-copy-numbers --input-allele-columns cn"))
}

include { MEDICC_SITKA } from './subworkflows/medicc_sitkasegs'

workflow {
    MEDICC_SITKA(signals, segments, medicc_args, allele_specific, cell_list, output_directory)
}
