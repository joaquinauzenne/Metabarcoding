#!/usr/bin/env python
# coding: utf-8

# Adapted from PIMENTA from Valerie van der Vorst et al., 2024

import os
import logging
from glob import glob
from io import StringIO

import click
import pandas as pd
from Bio import SeqIO

LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(filename)s:%(lineno)d]'
logging.basicConfig(format=LOG_FORMAT, level=logging.INFO)

@click.command()
@click.option('-s', '--sequencing-summary', required=True, help='Guppy basecaller sequencing_summary.txt output file')
@click.option('-b', '--barcoding-summary', required=True, help='Guppy barcoder barcoding_summary.txt output file')
@click.option('-i','--input-dir', required=True, help='Multiplexed FASTQ input directory')
@click.option('-o','--output-dir', required=True, help='Demultiplexed FASTQ output directory')

def demux(sequencing_summary, barcoding_summary, input_dir, output_dir):
    df_seq_summary = pd.read_table(sequencing_summary).set_index('read_id')
    logging.info(f'Loaded "{sequencing_summary}"')
    df_barcoding = pd.read_table(barcoding_summary).set_index('read_id')
    logging.info(f'Loaded "{barcoding_summary}"')

    n_runids = df_seq_summary.run_id.unique().size
    assert n_runids == 1, f'Expecting only one run_id for input dataset. Got {n_runids} ({df_seq_summary.run_id.value_counts()})'
    runid = df_seq_summary.run_id.unique()[0]

    try:
        os.makedirs(output_dir)
    except:
        logging.warning(f'Demux output directory "{output_dir} already exists!')

    fastqs = glob(os.path.join(input_dir, '*.fastq'))
    logging.info(f'Found {len(fastqs)} in "{input_dir}"')
    barcode_to_fh = None
    try:
        barcodes = df_barcoding.barcode_arrangement.unique()
        logging.info(f'Initializing file handles for {barcodes.size} barcodes')
        barcode_to_fh = {bc: open(os.path.join(output_dir, f'{bc}-{runid}.fq'), 'w') for bc in barcodes}
        count_demux = 0
        for fq in fastqs:
            rec_count = 0
            for rec in SeqIO.parse(fq, 'fastq'):
                rec_count += 1
                if rec.id in df_seq_summary.index and rec.id in df_barcoding.index:
                    if df_seq_summary.loc[rec.id, 'passes_filtering']:
                        bc = df_barcoding.loc[rec.id, 'barcode_arrangement']
                        sio = StringIO()
                        SeqIO.write([rec], sio, 'fastq')
                        barcode_to_fh[bc].write(sio.getvalue())
                        count_demux += 1
                else:
                    logging.warning(f'Read {rec.id} not in seq summary or barcoding results table! {rec.description}')
                if rec_count % 10000 == 0:
                    logging.info(f'Parsed {rec_count} of "{fq}". Demultiplexed {count_demux}')
    finally:
        if barcode_to_fh is not None and isinstance(barcode_to_fh, dict):
            logging.info(f'Closing file handles for {barcodes.size} barcodes')
            for fh in barcode_to_fh.values():
                fh.close()

if __name__ == '__main__':
    demux()