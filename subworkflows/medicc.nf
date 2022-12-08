
params.medicc_args = """-j 400 --input-type t --verbose --plot none --no-plot-tree \
--chromosomes-bed /juno/work/shah/isabl_software/dependencies/medicc2/medicc/objects/hg19_chromosome_arms.bed \
--regions-bed /juno/work/shah/users/myersm2/misseg/sitka-medicc-reconstruct/Davoli_2013_TSG_OG_genes_hg37.bed"""

process START_MEDICC_PARALLEL {
    input:
        tuple val(id), path(medicc_input), val(medicc_args)

    output:
        tuple val(id), path('task_*.pickle'), emit: tasks
        tuple val(id), path('tasks.pickle'), emit: task_idxs
        tuple val(id), path('sample_labels.pickle'), emit: sample_labels

    script:
    """
    medicc2 -j 400 --input-type t --verbose --plot none --no-plot-tree ${medicc_args} --start-external-parallel --task-dir ./ ${medicc_input} dummy_path
    """
}


process RUN_MEDICC_TASK {
    input:
        tuple val(id), path(t)
    
    output:
        tuple val(id), path("${t.baseName}_result.pickle"), emit: result

    script:
    """
    #!/usr/bin/env python

    import pickle

    with open('${t}', 'rb') as f:
        task = pickle.load(f)

    result = task[0](*task[1], **task[2])

    with open('${t.baseName}_result.pickle', 'wb') as f:
        pickle.dump(result, f)
    """
}


process FINISH_MEDICC_PARALLEL {
    input:
        tuple val(id), path(result), path(task_idxs), path(sample_labels), path(medicc_input), val(medicc_args), val(output_directory)

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
    medicc2 -j 400 --input-type t --verbose --plot none --no-plot-tree --chromosomes-bed \$CHROMOSOMES_BED --regions-bed \$REGIONS_BED ${medicc_args} --finish-external-parallel --task-dir ./ ${medicc_input} ./
    """
}


workflow MEDICC {
    take:
        medicc_input
        medicc_args
        output_directory
    main:
        START_MEDICC_PARALLEL(
            medicc_input
            .join(medicc_args))
        RUN_MEDICC_TASK(
            START_MEDICC_PARALLEL.out.tasks.transpose())
        FINISH_MEDICC_PARALLEL(
            RUN_MEDICC_TASK.out.result.groupTuple()
            .join(START_MEDICC_PARALLEL.out.task_idxs)
            .join(START_MEDICC_PARALLEL.out.sample_labels)
            .join(medicc_input)
            .join(medicc_args)
            .join(output_directory))
    emit:
        copynumber_events_df = FINISH_MEDICC_PARALLEL.out.copynumber_events_df
        events_overlap = FINISH_MEDICC_PARALLEL.out.events_overlap
        final_cn_profiles = FINISH_MEDICC_PARALLEL.out.final_cn_profiles
        final_tree_newick = FINISH_MEDICC_PARALLEL.out.final_tree_newick
        final_tree_xml = FINISH_MEDICC_PARALLEL.out.final_tree_xml
        pairwise_distances = FINISH_MEDICC_PARALLEL.out.pairwise_distances
        summary = FINISH_MEDICC_PARALLEL.out.summary
}

