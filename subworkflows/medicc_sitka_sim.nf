#!/usr/bin/env nextflow
params.medicc_args = """-j 400 --input-type t --verbose --plot none --no-plot-tree \
--chromosomes-bed /juno/work/shah/isabl_software/dependencies/medicc2/medicc/objects/hg19_chromosome_arms.bed \
--regions-bed /juno/work/shah/users/myersm2/misseg/sitka-medicc-reconstruct/Davoli_2013_TSG_OG_genes_hg37.bed"""

process CONVERT_SITKA_TREE {
    input:
        tuple val(id), path(tree), path(medicc_input), val(output_directory), val(cell_list)

    output:
        tuple val(id), path("${id}_tree_converted.newick"), emit: tree_converted
        tuple val(id), path("${id}_input_filtered.tsv"), emit: input_filtered

    publishDir "${output_directory}", mode: 'copy', overwrite: true

    script:
    """
    reformat_sitka_tree_sim.py ${tree} ${medicc_input} ${id}_tree_converted.newick ${id}_input_filtered.tsv --cell_list ${cell_list}
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
    medicc2 ${params.medicc_args} ${medicc_args} --tree ${tree} --events ${medicc_input} ./ --prefix ${id}
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

workflow MEDICC_SITKA_SIM {
    take:
        tree
        medicc_input
        medicc_args
        allele_specific
        cell_list
        output_directory
    main:
        CONVERT_SITKA_TREE(tree
            .join(medicc_input)
            .join(output_directory)
            .join(cell_list))
        RUN_MEDICC_WITH_TREE(CONVERT_SITKA_TREE.out.input_filtered
            .join(CONVERT_SITKA_TREE.out.tree_converted)
            .join(medicc_args)
            .join(output_directory))
        PLOT_MEDICC_RESULTS(medicc_input
            .join(RUN_MEDICC_WITH_TREE.out.final_cn_profiles)
            .join(RUN_MEDICC_WITH_TREE.out.final_tree_newick)
            .join(allele_specific)
            .join(output_directory))
    emit:
        tree_converted = CONVERT_SITKA_TREE.out.tree_converted
        copynumber_events_df = RUN_MEDICC_WITH_TREE.out.copynumber_events_df
        events_overlap = RUN_MEDICC_WITH_TREE.out.events_overlap
        final_cn_profiles = RUN_MEDICC_WITH_TREE.out.final_cn_profiles
        final_tree_newick = RUN_MEDICC_WITH_TREE.out.final_tree_newick
        final_tree_xml = RUN_MEDICC_WITH_TREE.out.final_tree_xml
        pairwise_distances = RUN_MEDICC_WITH_TREE.out.pairwise_distances
        summary = RUN_MEDICC_WITH_TREE.out.summary
        tree_cn_figure = PLOT_MEDICC_RESULTS.out.tree_cn_figure
}
