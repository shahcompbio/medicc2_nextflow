process RUN_MEDICC {
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
    CHROMOSOMES_BED=/juno/work/shah/isabl_software/dependencies/medicc2/medicc/objects/hg19_chromosome_arms.bed
    REGIONS_BED=/juno/work/shah/users/myersm2/misseg/sitka-medicc-reconstruct/Davoli_2013_TSG_OG_genes_hg37.bed 
    medicc2 -j 32 --input-type t --events --verbose --plot none --chromosomes-bed \$CHROMOSOMES_BED --regions-bed \$REGIONS_BED ${medicc_args} ${medicc_input} ./
    """
}

workflow MEDICC {
    take:
        medicc_input
        medicc_args
        output_directory
    main:
        RUN_MEDICC(
            medicc_input
            .join(medicc_args)
            .join(output_directory))
    emit:
        copynumber_events_df = RUN_MEDICC.out.copynumber_events_df
        events_overlap = RUN_MEDICC.out.events_overlap
        final_cn_profiles = RUN_MEDICC.out.final_cn_profiles
        final_tree_newick = RUN_MEDICC.out.final_tree_newick
        final_tree_xml = RUN_MEDICC.out.final_tree_xml
        pairwise_distances = RUN_MEDICC.out.pairwise_distances
        summary = RUN_MEDICC.out.summary
}

