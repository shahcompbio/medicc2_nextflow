
params.medicc_args = """-j 8 --input-type t --verbose --plot none --no-plot-tree --events \
--chromosomes-bed /rtsess01/compute/juno/shah/users/myersm2/repos/medicc2/medicc/objects/hg19_chromosome_arms.bed \
--regions-bed /rtsess01/compute/juno/shah/users/myersm2/reference/Davoli_2013_TSG_OG_genes_hg37.bed"""

process RUN_MEDICC_SERIAL {
    conda '/juno/home/myersm2/miniconda3/envs/medicc_force_clonal'

    input:
        tuple val(id), path(medicc_input), val(medicc_args), val(output_directory)

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
    medicc2 ${params.medicc_args} ${medicc_args} ${medicc_input} ./
    """

}


workflow MEDICC {
    take:
        medicc_input
        medicc_args
        output_directory
    main:
        RUN_MEDICC_SERIAL(
            medicc_input
            .join(medicc_args)
            .join(output_directory))
    emit:
        copynumber_events_df = RUN_MEDICC_SERIAL.out.copynumber_events_df
        events_overlap = RUN_MEDICC_SERIAL.out.events_overlap
        final_cn_profiles = RUN_MEDICC_SERIAL.out.final_cn_profiles
        final_tree_newick = RUN_MEDICC_SERIAL.out.final_tree_newick
        final_tree_xml = RUN_MEDICC_SERIAL.out.final_tree_xml
        pairwise_distances = RUN_MEDICC_SERIAL.out.pairwise_distances
        summary = RUN_MEDICC_SERIAL.out.summary
}

