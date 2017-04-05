from Bio import Entrez
import click
from progressbar import ProgressBar, UnknownLength
from zoo.load import write_record, accessions_to_fmt, read_accessions


# @click.option(
#     '--out', required=True,
#     help='The directory to write downloaded files to.')
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
    help='Output format: genbank | fasta | json')
@click.option(
    '--source', required=False, default='ncbi',
    help='The data source.')
@click.option(
    '--ids', required=True,
    help='A list of IDs, e.g. GenBank accession IDs.')
@click.argument('out', type=click.Path())
@click.command()
def load(ids, out, batch, email, db, fmt, source):
    '''Download genbank entries from NCBI and optionally dump to JSON

    This function borrows heavily from https://www.biostars.org/p/66921/,
    see full script at zoo/zoo/scripts/download_ncbi.py

    Help:

    out .. A file name (json) or directory (genbank, fasta).

    Example:

    zoo load --source ncbi --fmt json --email '' \
    --ids data/rna_virome_shi2016/rna_virome_shi2016.txt \
    result.json

    \b
    Loading data from NCBI.
    Batch size: 100
    Starting download:
    \ 10 Elapsed Time: 0:00:00
    Done.

    Notes:

    The function can fail due to a variety of connection problems with NCBI.
    If this is the case, get a coffee, and simply try again.

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

    For JSON as output format, data from each record is appended to the
    specified "out" file. If it exists already, data will get appended,
    modifying the existing file, which might be not what you want.
    '''
    if source == 'ncbi':
        print('Loading data from NCBI.')
        RETMAX = 10**9
        accessions = read_accessions(ids)
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
            write_record(out, acc, record, fmt)
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
