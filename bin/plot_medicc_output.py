#!/usr/bin/env python

import logging
import sys
import scgenome.loaders.annotation
import scgenome.loaders.hmmcopy
import scgenome.utils
import Bio.Phylo
import pandas as pd
import click
import matplotlib.pyplot as plt


def load_tree(newick_filename):
    tree = Bio.Phylo.read(newick_filename, 'newick')
    scgenome.tl.prune_leaves(tree, lambda a: a.name == 'diploid')
    return tree


medicc_input_dtype = {
    'sample_id': 'category',
    'chrom': 'category',
    'start': int,
    'end': int,
    'chr': 'category',
    'original_sample_id': 'category',
    'original_library_id': 'category',
    'cn_a': int,
    'cn_b': int,
}


def plot_tree_cn(medicc_input_filename, tree_filename, cn_profiles_filename, tree_cn_figure, allele_specific):
    medicc_input = pd.read_csv(medicc_input_filename, sep='\t', dtype=medicc_input_dtype)

    cell_info = medicc_input.rename(columns={
        'sample_id': 'cell_id',
        'original_sample_id': 'sample_id',
        'original_library_id': 'library_id',
    })[['cell_id', 'sample_id', 'library_id']].drop_duplicates()

    adata = scgenome.pp.load_cn.read_medicc2_cn(
        cn_profiles_filename, allele_specific=allele_specific)

    tree = load_tree(tree_filename)

    adata.obs['sample_id'] = cell_info.set_index('cell_id')['sample_id']

    tree, adata = scgenome.tl.align_cn_tree(tree, adata)

    scgenome.pl.plot_tree_cn(
        tree, adata,
        obs_annotation='sample_id', obs_cmap=plt.get_cmap("tab10"))

    plt.gcf().savefig(tree_cn_figure, bbox_inches='tight')


@click.command()
@click.argument('medicc_input_filename')
@click.argument('cn_profiles_filename')
@click.argument('tree_filename')
@click.argument('tree_cn_figure')
@click.option('--allele_specific', is_flag=True)
def medicc_plots(
        medicc_input_filename,
        cn_profiles_filename,
        tree_filename,
        tree_cn_figure,
        allele_specific=False,
    ):

    plot_tree_cn(medicc_input_filename, tree_filename, cn_profiles_filename, tree_cn_figure, allele_specific=allele_specific)


if __name__ == "__main__":
    logging.basicConfig(stream=sys.stderr, level=logging.INFO)
    medicc_plots()
