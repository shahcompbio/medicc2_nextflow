#!/usr/bin/env python

import logging
import sys
import Bio.Phylo
import pandas as pd
import click
import matplotlib.pyplot as plt
import anndata as ad
import numpy as np
import pyranges as pr

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

def dataframe_to_pyranges(data):
    data = pr.PyRanges(data.rename(columns={
        'chr': 'Chromosome',
        'start': 'Start',
        'end': 'End',
    }))

    return data

def pyranges_to_dataframe(data):
    data = data.as_df().rename(columns={
        'Chromosome': 'chr',
        'Start': 'start',
        'End': 'end',
    })

    return data

def segments2bins(seg_adata, binsize = int(5e5)):
    output_layers = seg_adata.layers
    
    segments = seg_adata.var
    binsize = int(5e5)

    bin_rows = []
    for _, r in seg_adata.var.iterrows():
        start = r.start
        while start < r.end:
            bin_rows.append([f'{r.chr}:{start}-{start + binsize - 1}', r.chr, start, start + binsize - 1])
            start += binsize
    bins = pd.DataFrame(bin_rows, columns = ['bin', 'chr', 'start', 'end']).set_index('bin').sort_values(by = ['chr', 'start'])
    
    # Match segments to bins
    bins['bin_index'] = np.arange(len(bins))
    segments['segment_index'] = np.arange(len(segments))

    pyr_bins = dataframe_to_pyranges(bins)
    pyr_segments = dataframe_to_pyranges(segments)

    intersect_1 = pyr_segments.intersect(pyr_bins)
    intersect_2 = pyr_bins.intersect(pyr_segments)

    intersect = pd.merge(
        pyranges_to_dataframe(intersect_1),
        pyranges_to_dataframe(intersect_2))

    bins = bins.merge(intersect, on = ['chr', 'start', 'end', 'bin_index'])
    
    # Resolve outputs: replicate out CN signal to original bin shape
    out_adata = ad.AnnData(np.zeros((seg_adata.shape[0], len(bins))),
                          obs = seg_adata.obs.copy(),
                          var = bins.copy())
        
    layers = {c:np.zeros(out_adata.shape).T for c in output_layers}
    newX = np.zeros(out_adata.shape).T 

    for seg_id, df in intersect.groupby('segment_index'):
        # select all bins that match this segment
        bin_idx = df.bin_index.to_numpy()
    
        for layer_id, arr in layers.items():
            arr[bin_idx] = seg_adata.layers[layer_id][:, seg_id].flatten()
        newX[bin_idx] =  seg_adata.X[:, seg_id].flatten()
    
    for layer_id, arr in layers.items():
        out_adata.layers[layer_id] = arr.T
    out_adata.X = newX.T
    
    return out_adata

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

    if not np.all(adata.var.end - adata.var.start == adata.var.iloc[0].end - adata.var.iloc[0].start):
        # Bins/segments have different lengths
        adata = segments2bins(adata)

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
