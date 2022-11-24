#!/usr/bin/env python

import Bio
import click
import re
import itertools
import random
import sys
import pandas as pd
from Bio import Phylo

sys.setrecursionlimit(1500)


def prune_leaves(clade, f):
    new_clades = []
    for c in clade.clades:
        c = prune_leaves(c, f)
        if not c.is_terminal() or not f(c):
            new_clades.append(c)
    clade.clades = new_clades
    return clade


def prune_leaves2(tree, f):
    while tree.count_terminals() > 0:
        # Find all terminals to be removed
        to_remove = []
        for n in tree.get_terminals():
            if f(n):
                to_remove.append(n.name)
        # Recreate clades without filtering terminals
        for n in tree.find_clades():
            new_clades = []
            for n2 in n.clades:
                if not n2.is_terminal() or not n2.name in to_remove:
                    new_clades.append(n2)
            n.clades = new_clades
        if len(to_remove) == 0:
            break
    return tree


def merge_single_child_clades(clade):
    if len(clade.clades) == 0:
        return clade
    elif len(clade.clades) == 1:
        new_clade = merge_single_child_clades(clade.clades[0])
        new_clade.branch_length += clade.branch_length
        return new_clade
    else:
        clade.clades = [merge_single_child_clades(c) for c in clade.clades]
        return clade


def merge_single_child_clades2(tree):
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

    return tree


def split_polytomies(clade):
    if len(clade.clades) > 2:
        new_clades = []
        for c1, c2 in  list(itertools.zip_longest(clade.clades[::2], clade.clades[1::2])):
            if c2 is not None:
                new_clade = Bio.Phylo.BaseTree.Clade(branch_length=0.)
                new_clade.clades = [c1, c2]
                new_clades.append(new_clade)
            else:
                new_clades.append(c1)
        clade.clades = new_clades
        return split_polytomies(clade)
    elif len(clade.clades) <= 2:
        clade.clades = [split_polytomies(c) for c in clade.clades]
        return clade


def split_polytomies2(tree):
    while tree.count_terminals() > 0:
        altered = False
        for n in tree.find_clades():
            new_clades = []
            if len(n.clades) > 2:
                for c1, c2 in  list(itertools.zip_longest(n.clades[::2], n.clades[1::2])):
                    if c2 is not None:
                        new_clade = Bio.Phylo.BaseTree.Clade(branch_length=0.)
                        new_clade.clades = [c1, c2]
                        new_clades.append(new_clade)
                    else:
                        new_clades.append(c1)
                n.clades = new_clades
                altered = True
                break
        if not altered:
            break
    
    return tree


def test_merge_single_child_clades(tree_filename):
    tree1 = Phylo.read(tree_filename, 'newick')
    tree2 = Phylo.read(tree_filename, 'newick')

    # Add branch lengths
    for n in tree1.find_clades():
        n.branch_length = 1.

    # Add branch lengths
    for n in tree2.find_clades():
        n.branch_length = 1.

    tree1 = Bio.Phylo.BaseTree.Tree.from_clade(merge_single_child_clades(tree1.clade))
    tree1.rooted = False
    tree1.weight = 1.0

    tree2 = merge_single_child_clades2(tree2)
    tree2.rooted = False
    tree2.weight = 1.0

    assert str(tree1) == str(tree2)


def test_split_polytomies(tree_filename):
    tree1 = Phylo.read(tree_filename, 'newick')
    tree2 = Phylo.read(tree_filename, 'newick')

    tree1 = Bio.Phylo.BaseTree.Tree.from_clade(split_polytomies(tree1.clade))
    tree1.rooted = False
    tree1.weight = 1.0

    tree2 = split_polytomies2(tree2)
    tree2.rooted = False
    tree2.weight = 1.0

    assert str(tree1) == str(tree2)


def test_prune_leaves(tree_filename):
    tree1 = Phylo.read(tree_filename, 'newick')
    tree2 = Phylo.read(tree_filename, 'newick')

    selected = [n.name for n in tree1.get_terminals()]
    random.shuffle(selected)
    selected = selected[:3]

    tree1 = Bio.Phylo.BaseTree.Tree.from_clade(prune_leaves(tree1.clade, lambda n: n.name not in selected))
    tree1.rooted = False
    tree1.weight = 1.0

    tree2 = prune_leaves2(tree2, lambda n: n.name not in selected)
    tree2.rooted = False
    tree2.weight = 1.0

    assert str(tree1) == str(tree2)


def test_tree_algs(tree_filename):
    random.seed(1)
    test_merge_single_child_clades(tree_filename)
    test_split_polytomies(tree_filename)
    test_prune_leaves(tree_filename)


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
    tree = Bio.Phylo.BaseTree.Tree.from_clade(prune_leaves(tree.clade, lambda n: n.name not in intersection_cells))

    # Remove chains
    tree = Bio.Phylo.BaseTree.Tree.from_clade(merge_single_child_clades(tree.clade))

    # Binarize the tree
    tree = Bio.Phylo.BaseTree.Tree.from_clade(split_polytomies(tree.clade))

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
    test_tree_algs(tree_filename)
    prepare_inputs(tree_filename, table_filename, reformatted_tree_filename, reformatted_table_filename)


if __name__ == "__main__":
    reformat_sitka_tree_for_medicc()
