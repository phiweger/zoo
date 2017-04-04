import Bio.Entrez
import click
from progressbar import ProgressBar, UnknownLength
import os
from zoo.load import write_record, accessions_to_gb, read_accessions


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
@click.argument('file', type=click.Path())
@click.command()
def load(file, out, batch, email, db):
    '''Download genbank entries from NCBI and optionally dump to JSON

    This function borrows heavily from https://www.biostars.org/p/66921/.

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
    '''
    print('Loading data from NCBI.')

    RETMAX = 10**9

    accessions = read_accessions(file)
    op_dir = out
    if not os.path.exists(op_dir):
        os.makedirs(op_dir)
    dbase = db
    Bio.Entrez.email = email
    batchsize = batch

    print('Batch size:', batchsize)
    print('Starting download:')
    bar = ProgressBar(max_value=UnknownLength)
    counter = 0
    for acc, record in accessions_to_gb(
            accessions, dbase, batchsize, RETMAX):
        counter += 1
        write_record(op_dir, acc, record)
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
