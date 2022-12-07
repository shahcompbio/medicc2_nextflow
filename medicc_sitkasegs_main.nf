#!/usr/bin/env nextflow

params.id ?: { params.id = 'output' }()

params.signals_filename ?: { log.error "No signals data provided. Make sure you have used the '--signals_filename' option."; exit 1 }()
signals = Channel.of(tuple(params.id, file(params.signals_filename, checkIfExists: true)))

params.segments_filename ?: { log.error "No segments provided. Make sure you have used the '--segments_filename' option."; exit 1 }()
segments = Channel.of(tuple(params.id, file(params.segments_filename, checkIfExists: true)))

params.medicc_args = "-j 400 --input-type t --verbose --plot none --no-plot-tree"
medicc_args = Channel.of(tuple(params.id, params.medicc_args))

params.allele_specific = "--allele_specific"
allele_specific = Channel.of(tuple(params.id, params.allele_specific))

params.output_directory ?: { log.error "No copy number data provided. Make sure you have used the '--output_directory' option."; exit 1 }()
output_directory = Channel.of(tuple(params.id, params.output_directory))

include { MEDICC_SITKA } from './subworkflows/medicc_sitka'

workflow {
    MEDICC_SITKA(signals, segments, medicc_args, allele_specific, output_directory)
}
