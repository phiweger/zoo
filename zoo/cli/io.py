from Bio import Entrez
import click
from progressbar import ProgressBar, UnknownLength
import os
from zoo.load import write_record, accessions_to_fmt, read_accessions


@click.option(
    '--out', required=True,
    help='The directory to write downloaded files to.')
@click.option(
    '--batch', required=False, default=100,
    help='The number of accessions to process per request.')
@click.option(
    '--email', required=False, default='',
    help='An e-mail address.')
@click.option(
    '--db', required=False, default='nucleotide',
    help='"NCBI database ID.')
@click.option(
    '--fmt', required=False, default='genbank',
    help='"Output format: genbank | fasta | json')
@click.argument('file', type=click.Path())
@click.command()
def load(file, out, batch, email, db, fmt):
    '''Download genbank entries from NCBI and optionally dump to JSON

    This function borrows heavily from https://www.biostars.org/p/66921/,
    see full script at zoo/zoo/scripts/download_ncbi.py

    Help:

    file .. A file with accessions to download.

    Desired syntax:

    \b
    zoo load --ncbi accession.txt result.json

    Example:

    zoo load --out virometest6 --email '' rna_virome_shi2016.txt

    \b
    Loading data from NCBI.
    Batch size: 100
    Starting download:
    \ 10 Elapsed Time: 0:00:00
    Done.

    Notes:

    The function can fail with e.g.:
    urllib.error.URLError: <urlopen error [Errno 54] Connection reset by peer>
    probably due to some connection problem with NCBI. If this is the case,
    get a coffee, and simply try again.

    \b
    Loading data from NCBI.
    Batch size: 100
    Starting download:
    / 200 Elapsed Time: 0:00:08

    \b
    <urlopen error [Errno 60] Operation timed out>
    Reconnect, attempt 1/10.
    / 600 Elapsed Time: 0:02:29

    \b
    <urlopen error [Errno 54] Connection reset by peer>
    Reconnect, attempt 1/10.
    - 2109 Elapsed Time: 0:07:06
    Done.
    '''
    print('Loading data from NCBI.')

    RETMAX = 10**9

    accessions = read_accessions(file)
    if not os.path.exists(out):
        os.makedirs(out)
    dbase = db
    Entrez.email = email
    batchsize = batch

    print('Batch size:', batchsize)
    print('Starting download:')
    bar = ProgressBar(max_value=UnknownLength)
    counter = 0
    for acc, record in accessions_to_fmt(
            accessions, dbase, batchsize, RETMAX, fmt):
        counter += 1
        write_record(out, acc, record)
        bar.update(counter)
    print('\nDone.')


@click.command()
def dump():
    '''
    Example:

    \b
    zoo dump --db x --cell y --fasta dump.fa  # all metadata in header, | delim
    '''
    print('Dump.')
