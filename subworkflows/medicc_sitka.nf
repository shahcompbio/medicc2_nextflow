#!/usr/bin/env nextflow


process CONVERT_SITKA_TREE {
    input:
        tuple val(id), path(tree), path(signals), val(output_directory)

    output:
        tuple val(id), path("${id}_tree_converted.newick"), emit: tree_converted
        tuple val(id), path("${id}_signals_filtered.csv.gz"), emit: signals_filtered

    publishDir "${output_directory}", mode: 'copy', overwrite: true

    script:
    """
    reformat_sitka_tree.py ${tree} ${signals} ${id}_tree_converted.newick ${id}_signals_filtered.csv.gz
    """
}

process GENERATE_MEDICC_INPUT {
    input:
        tuple val(id), path(signals), val(allele_specific), val(output_directory)

    output:
        tuple val(id), path("${id}.tsv"), emit: medicc_input

    publishDir "${output_directory}", mode: 'copy', overwrite: true

    script:
    """
    create_medicc_input.py ${id}.tsv --signals_results ${signals} ${allele_specific}
    """
}

process RUN_MEDICC_WITH_TREE {
    input:
        tuple val(id), path(medicc_input), path(tree), val(medicc_args), val(output_directory)

    output:
        tuple val(id), path("${id}_copynumber_events_df.tsv"), emit: copynumber_events_df
        tuple val(id), path("${id}_events_overlap.tsv"), emit: events_overlap
        tuple val(id), path("${id}_final_cn_profiles.tsv"), emit: final_cn_profiles
        tuple val(id), path("${id}_final_tree.new"), emit: final_tree_newick
        tuple val(id), path("${id}_final_tree.xml"), emit: final_tree_xml
        tuple val(id), path("${id}_pairwise_distances.tsv"), emit: pairwise_distances
        tuple val(id), path("${id}_summary.tsv"), emit: summary

    publishDir "${output_directory}", mode: 'copy', overwrite: true

    script:
    """
    medicc2 -j 400 --input-type t --verbose --plot none --no-plot-tree ${medicc_args} --tree ${tree} ${medicc_input} ./
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
        tree
        signals
        medicc_args
        allele_specific
        output_directory
    main:
        CONVERT_SITKA_TREE(tree
            .join(signals)
            .join(output_directory))
        GENERATE_MEDICC_INPUT(CONVERT_SITKA_TREE.out.signals_filtered
            .join(allele_specific)
            .join(output_directory))
        RUN_MEDICC_WITH_TREE(GENERATE_MEDICC_INPUT.out.medicc_input
            .join(CONVERT_SITKA_TREE.out.tree_converted)
            .join(medicc_args)
            .join(output_directory))
        PLOT_MEDICC_RESULTS(GENERATE_MEDICC_INPUT.out.medicc_input
            .join(RUN_MEDICC_WITH_TREE.out.final_cn_profiles)
            .join(RUN_MEDICC_WITH_TREE.out.final_tree_newick)
            .join(allele_specific)
            .join(output_directory))
    emit:
        tree_converted = CONVERT_SITKA_TREE.out.tree_converted
        medicc_input = GENERATE_MEDICC_INPUT.out.medicc_input
        copynumber_events_df = RUN_MEDICC_WITH_TREE.out.copynumber_events_df
        events_overlap = RUN_MEDICC_WITH_TREE.out.events_overlap
        final_cn_profiles = RUN_MEDICC_WITH_TREE.out.final_cn_profiles
        final_tree_newick = RUN_MEDICC_WITH_TREE.out.final_tree_newick
        final_tree_xml = RUN_MEDICC_WITH_TREE.out.final_tree_xml
        pairwise_distances = RUN_MEDICC_WITH_TREE.out.pairwise_distances
        summary = RUN_MEDICC_WITH_TREE.out.summary
        tree_cn_figure = PLOT_MEDICC_RESULTS.out.tree_cn_figure
}
