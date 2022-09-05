#!/usr/bin/env nextflow

params.signals_results ?: { log.error "No signals data provided. Make sure you have used the '--signals_results' option."; exit 1 }()
signals_results = file(params.signals_results, checkIfExists: true)

params.annotation_metrics_csv ?: { log.error "No annotations data provided. Make sure you have used the '--annotation_metrics_csv' option."; exit 1 }()

annotation_metrics_csv = Channel.fromPath(params.annotation_metrics_csv) | splitCsv(header:true) | map { row -> file(row.annotation_metrics, checkIfExists: true) } | toList()
annotation_metrics_csv_yaml = Channel.fromPath(params.annotation_metrics_csv) | splitCsv(header:true) | map { row -> file(row.annotation_metrics + '.yaml', checkIfExists: true) } | toList()

params.output_directory ?: { log.error "No copy number data provided. Make sure you have used the '--output_directory' option."; exit 1 }()

params.results_basename ?: { params.results_basename = 'output' }()

params.medicc_args = "-j 400 --input-type t --verbose --plot none --no-plot-tree"

include { MEDICC } from './subworkflows/medicc'

process GENERATE_MEDICC_INPUT {
    input:
        path signals_results
        path annotation_metrics_csv, stageAs: 'annotation_metrics*.csv.gz'
        path annotation_metrics_csv_yaml, stageAs: 'annotation_metrics*.csv.gz.yaml'

    output:
        path 'medicc_input.tsv', emit: medicc_input

    script:
    """
    create_medicc_input.py medicc_input.tsv --signals_results ${signals_results} --annotation_metrics ${annotation_metrics_csv.join(' --annotation_metrics ')}
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
    plot_medicc_output.py ${medicc_input} ${cn_profiles} ${tree} ${params.results_basename}_tree_heatmap.png --allele_specific
    """
}

workflow {
    GENERATE_MEDICC_INPUT(signals_results, annotation_metrics_csv, annotation_metrics_csv_yaml)
    MEDICC(GENERATE_MEDICC_INPUT.out.medicc_input)
    PLOT_MEDICC_RESULTS(GENERATE_MEDICC_INPUT.out.medicc_input, MEDICC.out.final_cn_profiles, MEDICC.out.final_tree_newick)
    PLOT_MEDICC_RESULTS.out.tree_cn_figure.subscribe { it.copyTo(params.output_directory) }
}

