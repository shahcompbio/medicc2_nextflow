#!/usr/bin/env python

import logging
import sys
import scgenome.loaders.annotation
import scgenome.loaders.hmmcopy
import scgenome.utils
import pandas as pd
import click


@click.command()
@click.argument('output_filename')
@click.option('--hmmcopy_reads', multiple=True)
@click.option('--signals_results', multiple=True)
@click.option('--annotation_metrics', multiple=True)
@click.option('--allele_specific', is_flag=True)
def create_medicc2_input(
        output_filename,
        hmmcopy_reads,
        signals_results,
        annotation_metrics,
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
