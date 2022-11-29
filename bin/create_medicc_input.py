#!/usr/bin/env python

import logging
import sys
import scgenome.loaders.annotation
import scgenome.loaders.hmmcopy
import scgenome.utils
import pandas as pd
import pyranges as pr
import numpy as np
import click



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


def resegment(cn_data, segments, cn_cols):
    def create_segments(df):
        positions = np.unique(np.concatenate([(df['start'] - 1).values, df['end'].values]))
        return pd.DataFrame({'start': positions[:-1:] + 1, 'end': positions[1::]})

    # Consolodate segments
    segments = segments.groupby(['chr'], observed=True).apply(create_segments).reset_index()[['chr', 'start', 'end']].sort_values(['chr', 'start'])

    bins = cn_data[['chr', 'start', 'end']].drop_duplicates()

    segments['segment_idx'] = range(len(segments.index))
    bins['bin_idx'] = range(len(bins.index))

    pyr_bins = dataframe_to_pyranges(bins)
    pyr_segments = dataframe_to_pyranges(segments)

    intersect_1 = pyr_segments.intersect(pyr_bins)
    intersect_2 = pyr_bins.intersect(pyr_segments)

    intersect = pd.merge(
        pyranges_to_dataframe(intersect_1),
        pyranges_to_dataframe(intersect_2))

    cn_data = cn_data.merge(intersect, how='left')
    assert not cn_data['segment_idx'].isnull().any()

    segment_data = cn_data.groupby(['cell_id', 'segment_idx'], observed=True)[cn_cols].mean().round().astype(int).reset_index()
    segment_data = segment_data.merge(segments)

    return segment_data


@click.command()
@click.argument('output_filename')
@click.option('--hmmcopy_reads', multiple=True)
@click.option('--signals_results', multiple=True)
@click.option('--annotation_metrics', multiple=True)
@click.option('--segments_filename')
@click.option('--allele_specific', is_flag=True)
def create_medicc2_input(
        output_filename,
        hmmcopy_reads,
        signals_results,
        annotation_metrics,
        segments_filename,
        allele_specific,
    ):

    if len(hmmcopy_reads) > 0:
        assert len(signals_results) == 0
        assert len(annotation_metrics) > 0, 'medicc from raw hmmcopy will not work without filtering from annotation metrics'
        assert allele_specific is not True, 'can only do total copy number from hmmcopy input'

    elif len(signals_results) > 0:
        assert len(hmmcopy_reads) == 0

    else:
        raise ValueError('specify either signals or hmmcopy input')

    if allele_specific:
        cn_cols = ['cn_a', 'cn_b']
    
    else:
        cn_cols = ['cn']

    if len(hmmcopy_reads) > 0:
        cn_data = []
        for filename in hmmcopy_reads:
            cn_data.append(scgenome.loaders.hmmcopy.process_hmmcopy_data(
                filename, usecols=scgenome.loaders.hmmcopy.standard_hmmcopy_reads_cols))
        cn_data = scgenome.utils.concat_with_categories(cn_data, ignore_index=True)

        cn_data['cn'] = cn_data['state']

    else:
        signals_dtype = {
            'chr': 'category',
            'cell_id': 'category',
        }
        signals_cols = [
            'chr',
            'start',
            'end',
            'cell_id',
            'state',
            'Maj',
            'Min',
        ]
        cn_data = []
        for filename in signals_results:
            cn_data.append(pd.read_csv(filename, dtype=signals_dtype, usecols=signals_cols))
        cn_data = scgenome.utils.concat_with_categories(cn_data, ignore_index=True)

        cn_data['cn'] = cn_data['state']
        cn_data['cn_a'] = cn_data['Maj']
        cn_data['cn_b'] = cn_data['Min']

        # Remove chromosomes for which all bins are null in all cells (eg Y chromosome)
        non_null_chroms = (
            cn_data.set_index(['chr', 'start', 'end', 'cell_id'])[cn_cols].unstack()
                .isna().all(axis=1).groupby('chr').all().rename('null_chrom').reset_index().query('null_chrom == False'))
        cn_data = cn_data.merge(non_null_chroms[['chr']])

    # Resegment according to input segments
    if segments_filename is not None:
        segments = pd.read_csv(segments_filename, dtype={'chr': str})

        # Require the same set of chromosomes between segments and cn_data
        segments = segments[segments['chr'].isin(cn_data['chr'].unique())]
        cn_data = cn_data[cn_data['chr'].isin(segments['chr'].unique())]

        cn_data = resegment(cn_data, segments, cn_cols)

    else:
        raise ValueError('failed check')

    # Filter using annotation metrics if provided
    if len(annotation_metrics) > 0:
        metrics_data = []
        for filename in annotation_metrics:
            metrics_data.append(scgenome.loaders.annotation.process_annotation_file(filename))
        metrics_data = scgenome.utils.concat_with_categories(metrics_data, ignore_index=True)

        metrics_data = metrics_data.query('quality > 0.75 and not is_control')
        cn_data = cn_data.merge(
            metrics_data[['cell_id', 'library_id', 'sample_id']].rename(columns={
                'sample_id': 'original_sample_id',
                'library_id': 'original_library_id',
        }))

    # HACK: Parse cell id for sample and library
    else:
        cn_data['original_sample_id'] = cn_data['cell_id'].str.rsplit('-', expand=True, n=3)[0]
        cn_data['original_library_id'] = cn_data['cell_id'].str.rsplit('-', expand=True, n=3)[1]

    null_bins = (
        cn_data.set_index(['chr', 'start', 'end', 'cell_id'])[cn_cols].unstack()
            .isna().stack().any(axis=1).rename('null_bin').reset_index().query('null_bin == True'))

    if not null_bins.empty:
        raise ValueError(f'missing data for {null_bins.iloc[0].to_dict()}')

    medicc2_input = cn_data.rename(columns={'cell_id': 'sample_id'})

    autosomes = []
    for chr in cn_data['chr'].unique():
        try:
            autosomes.append(int(chr))
        except ValueError:
            if chr not in ['X', 'Y']:
                logging.warning(f'skipping chromosome {chr}')
            continue

    max_autosome = max(autosomes)
    chrom_X = max_autosome + 1
    chrom_Y = max_autosome + 2
    chromosomes_str = [str(a) for a in autosomes + ['X', 'Y']]

    medicc2_input = medicc2_input[medicc2_input['chr'].isin(chromosomes_str)]

    medicc2_input['chrom'] = medicc2_input['chr'].astype(str)
    medicc2_input.loc[medicc2_input['chr'] == 'X', 'chrom'] = chrom_X
    medicc2_input.loc[medicc2_input['chr'] == 'Y', 'chrom'] = chrom_Y
    medicc2_input['chrom'] = medicc2_input['chrom'].astype(int)

    cols = [
        'sample_id',
        'chrom',
        'start',
        'end',
        'chr',
        'original_sample_id',
        'original_library_id',
    ]

    medicc2_input[cols + cn_cols].to_csv(output_filename, sep='\t', index=False)


if __name__ == "__main__":
    logging.basicConfig(stream=sys.stderr, level=logging.INFO)
    create_medicc2_input()
