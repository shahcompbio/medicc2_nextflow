
params.medicc_args = """-j 400 --input-type t --verbose --plot none --no-plot-tree \
--chromosomes-bed /juno/work/shah/isabl_software/dependencies/medicc2/medicc/objects/hg19_chromosome_arms.bed \
--regions-bed /juno/work/shah/users/myersm2/misseg/sitka-medicc-reconstruct/Davoli_2013_TSG_OG_genes_hg37.bed"""

process START_MEDICC_PARALLEL {
    input:
        path "${params.results_basename}.tsv"

    output:
        path 'task_*.pickle', emit: tasks
        path 'tasks.pickle', emit: task_idxs
        path 'sample_labels.pickle', emit: sample_labels

    script:
    """
    medicc2 ${params.medicc_args} --start-external-parallel --task-dir ./ ${params.results_basename}.tsv dummy_path
    """
}

process RUN_MEDICC_TASK {
    input:
        path t
    
    output:
        path "${t.baseName}_result.pickle", emit: result

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
        path results
        path 'tasks.pickle'
        path 'sample_labels.pickle'
        path "${params.results_basename}.tsv"

    output:
        path "${params.results_basename}_copynumber_events_df.tsv", emit: copynumber_events_df
        path "${params.results_basename}_events_overlap.tsv", emit: events_overlap
        path "${params.results_basename}_final_cn_profiles.tsv", emit: final_cn_profiles
        path "${params.results_basename}_final_tree.new", emit: final_tree_newick
        path "${params.results_basename}_final_tree.xml", emit: final_tree_xml
        path "${params.results_basename}_pairwise_distances.tsv", emit: pairwise_distances
        path "${params.results_basename}_summary.tsv", emit: summary

    publishDir "${params.output_directory}", mode: 'copy'

    script:
    """
    medicc2 ${params.medicc_args} --finish-external-parallel --task-dir ./ ${params.results_basename}.tsv ./
    """
}


workflow MEDICC {
    take:
        copy_number
    main:
        START_MEDICC_PARALLEL(copy_number)
        RUN_MEDICC_TASK(START_MEDICC_PARALLEL.out.tasks.flatten())
        FINISH_MEDICC_PARALLEL(RUN_MEDICC_TASK.out.result.toList(), START_MEDICC_PARALLEL.out.task_idxs, START_MEDICC_PARALLEL.out.sample_labels, copy_number)
    emit:
        copynumber_events_df = FINISH_MEDICC_PARALLEL.out.copynumber_events_df
        events_overlap = FINISH_MEDICC_PARALLEL.out.events_overlap
        final_cn_profiles = FINISH_MEDICC_PARALLEL.out.final_cn_profiles
        final_tree_newick = FINISH_MEDICC_PARALLEL.out.final_tree_newick
        final_tree_xml = FINISH_MEDICC_PARALLEL.out.final_tree_xml
        pairwise_distances = FINISH_MEDICC_PARALLEL.out.pairwise_distances
        summary = FINISH_MEDICC_PARALLEL.out.summary
}

