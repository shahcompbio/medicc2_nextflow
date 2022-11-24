#!/usr/bin/env python

import Bio
import click
import re
import pandas as pd
from Bio import Phylo


def prune_leaves(tree, f):
    while tree.count_terminals() > 0:
        altered = False
        for a in tree.get_terminals():
            if f(a):
                tree.prune(a)
                altered = True
        if not altered:
            break


def prepare_inputs(tree_filename, table_filename, reformatted_tree_filename, reformatted_table_filename):
    # read in table and tree
    tree = Phylo.read(tree_filename, 'newick')
    table = pd.read_csv(table_filename)

    # Rename cell nodes
    for n in tree.find_clades():
        if n.name is not None:
            n.name = re.sub('^cell_', '', n.name)

    tree_cells = set([n.name for n in tree.get_terminals()])
    table_cells = set(table['cell_id'].unique())

    # restrict both to their intersection
    intersection_cells = tree_cells.intersection(table_cells)

    # Add branch lengths
    for n in tree.find_clades():
        n.branch_length = 1.

    table = table[table['cell_id'].isin(intersection_cells)]
    prune_leaves(tree, lambda n: n.name not in intersection_cells)

    # Merge single branch nodes with parent
    while tree.count_terminals() > 0:
        altered = False
        for n in tree.find_clades():
            new_clades = []
            for n2 in n.clades:
                if len(n2.clades) == 1:
                    n2_new = n2.clades[0]
                    n2_new.branch_length = n2.branch_length + n2.clades[0].branch_length
                    new_clades.append(n2_new)
                    altered = True
                else:
                    new_clades.append(n2)
            n.clades = new_clades
        if not altered:
            break

    # Merge single branch top clade if necessary
    if len(tree.clade.clades) == 1:
        new_clade = tree.clade.clades[0]
        new_clade.branch_length = new_clade.branch_length + tree.clade.branch_length
        tree = Bio.Phylo.BaseTree.Tree.from_clade(new_clade)

    # Split polytomies
    internal_idx = 0
    while tree.count_terminals() > 0:
        altered = False
        for n in tree.find_clades():
            if len(n.clades) > 2:
                new_subclade = Bio.Phylo.BaseTree.Clade(name=f'internal_extra_{internal_idx}', branch_length=0.)
                internal_idx += 1
                new_subclade.clades = n.clades[:2]
                n.clades = [new_subclade] + n.clades[2:]
                altered = True
                break
        if not altered:
            break

    # Relabel internal
    internal_idx = 0
    for n in tree.find_clades():
        assert len(n.clades) in [0, 2]
        if len(n.clades) == 2:
            n.name = f'internal_{internal_idx}'
            internal_idx += 1

    # Add dummy_normal node
    tree = Bio.Phylo.BaseTree.Tree.from_clade(
        Bio.Phylo.BaseTree.Clade(
            clades=[tree.clade, Bio.Phylo.BaseTree.Clade(name='diploid')]
        )
    )

    Bio.Phylo.write(tree, reformatted_tree_filename, 'newick')

    table.to_csv(reformatted_table_filename, index=False)


@click.command()
@click.argument('tree_filename')
@click.argument('table_filename')
@click.argument('reformatted_tree_filename')
@click.argument('reformatted_table_filename')
def reformat_sitka_tree_for_medicc(
        tree_filename,
        table_filename,
        reformatted_tree_filename,
        reformatted_table_filename
):
    prepare_inputs(tree_filename, table_filename, reformatted_tree_filename, reformatted_table_filename)


if __name__ == "__main__":
    reformat_sitka_tree_for_medicc()
