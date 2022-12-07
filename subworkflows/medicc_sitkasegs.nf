#!/usr/bin/env nextflow


include { MEDICC } from './medicc'


process GENERATE_MEDICC_INPUT {
    input:
        tuple val(id), path(signals), path(segments), val(allele_specific), val(output_directory)

    output:
        tuple val(id), path("${id}.tsv"), emit: medicc_input

    publishDir "${output_directory}", mode: 'copy', overwrite: true

    script:
    """
    create_medicc_input.py ${id}.tsv --signals_results ${signals} --segments_filename ${segments} ${allele_specific}
    """
}

process PLOT_MEDICC_RESULTS {
    input:
        tuple val(id), path(medicc_input), path(cn_profiles), path(tree), val(allele_specific), val(output_directory)

    output:
        tuple val(id), path("${id}_tree_heatmap.png"), emit: tree_cn_figure

    publishDir "${output_directory}", mode: 'copy', overwrite: true

    script:
    """
    plot_medicc_output.py ${medicc_input} ${cn_profiles} ${tree} ${id}_tree_heatmap.png ${allele_specific}
    """
}

workflow MEDICC_SITKA {
    take:
        signals
        segments
        medicc_args
        allele_specific
        output_directory
    main:
        GENERATE_MEDICC_INPUT(signals
            .join(segments)
            .join(allele_specific)
            .join(output_directory))
        MEDICC(
            GENERATE_MEDICC_INPUT.out.medicc_input,
            medicc_args,
            output_directory)
        PLOT_MEDICC_RESULTS(GENERATE_MEDICC_INPUT.out.medicc_input
            .join(MEDICC.out.final_cn_profiles)
            .join(MEDICC.out.final_tree_newick)
            .join(allele_specific)
            .join(output_directory))
    emit:
        medicc_input = GENERATE_MEDICC_INPUT.out.medicc_input
        copynumber_events_df = MEDICC.out.copynumber_events_df
        events_overlap = MEDICC.out.events_overlap
        final_cn_profiles = MEDICC.out.final_cn_profiles
        final_tree_newick = MEDICC.out.final_tree_newick
        final_tree_xml = MEDICC.out.final_tree_xml
        pairwise_distances = MEDICC.out.pairwise_distances
        summary = MEDICC.out.summary
        tree_cn_figure = PLOT_MEDICC_RESULTS.out.tree_cn_figure
}
