#!/usr/bin/env nextflow

params.hmmcopy_results_csv ?: { log.error "No annotations data provided. Make sure you have used the '--hmmcopy_results_csv' option."; exit 1 }()

hmmcopy_results_csv = Channel.fromPath(params.hmmcopy_results_csv) | splitCsv(header:true) | map { row -> file(row.hmmcopy_results, checkIfExists: true) } | toList()
hmmcopy_results_csv_yaml = Channel.fromPath(params.hmmcopy_results_csv) | splitCsv(header:true) | map { row -> file(row.hmmcopy_results + '.yaml', checkIfExists: true) } | toList()

params.annotation_metrics_csv ?: { log.error "No annotations data provided. Make sure you have used the '--annotation_metrics_csv' option."; exit 1 }()

annotation_metrics_csv = Channel.fromPath(params.annotation_metrics_csv) | splitCsv(header:true) | map { row -> file(row.annotation_metrics, checkIfExists: true) } | toList()
annotation_metrics_csv_yaml = Channel.fromPath(params.annotation_metrics_csv) | splitCsv(header:true) | map { row -> file(row.annotation_metrics + '.yaml', checkIfExists: true) } | toList()

params.output_directory ?: { log.error "No copy number data provided. Make sure you have used the '--output_directory' option."; exit 1 }()

params.results_basename ?: { params.results_basename = 'output' }()

params.medicc_args = "-j 400 --input-type t --total-copy-numbers --input-allele-columns cn --verbose --plot none --no-plot-tree"

include { MEDICC } from './subworkflows/medicc'

process GENERATE_MEDICC_INPUT {
    input:
        path hmmcopy_results_csv, stageAs: 'hmmcopy_results*.csv.gz'
        path hmmcopy_results_csv_yaml, stageAs: 'hmmcopy_results*.csv.gz.yaml'
        path annotation_metrics_csv, stageAs: 'annotation_metrics*.csv.gz'
        path annotation_metrics_csv_yaml, stageAs: 'annotation_metrics*.csv.gz.yaml'

    output:
        path "${params.results_basename}.tsv", emit: medicc_input

    script:
    """
    create_medicc_input.py ${params.results_basename}.tsv --hmmcopy_reads ${hmmcopy_results_csv.join(' --hmmcopy_reads ')} --annotation_metrics ${annotation_metrics_csv.join(' --annotation_metrics ')}
    """
}

process PLOT_MEDICC_RESULTS {
    input:
        path medicc_input
        path cn_profiles
        path tree

    output:
        path "${params.results_basename}_tree_heatmap.png", emit: tree_cn_figure

    script:
    """
    plot_medicc_output.py ${medicc_input} ${cn_profiles} ${tree} ${params.results_basename}_tree_heatmap.png
    """
}

workflow {
    GENERATE_MEDICC_INPUT(hmmcopy_results_csv, hmmcopy_results_csv_yaml, annotation_metrics_csv, annotation_metrics_csv_yaml)
    MEDICC(GENERATE_MEDICC_INPUT.out.medicc_input)
    PLOT_MEDICC_RESULTS(GENERATE_MEDICC_INPUT.out.medicc_input, MEDICC.out.final_cn_profiles, MEDICC.out.final_tree_newick)

    GENERATE_MEDICC_INPUT.out.medicc_input.subscribe { it.copyTo(params.output_directory) }
    PLOT_MEDICC_RESULTS.out.tree_cn_figure.subscribe { it.copyTo(params.output_directory) }
}

