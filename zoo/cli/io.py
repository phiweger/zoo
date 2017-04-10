from Bio import Entrez
import click
import json
from zoo.load import write_record, accessions_to_fmt, read_accessions
from zoo.utils import eprint, deep_get
from zoo.utils_mongodb import env_switch, eval_query
import sys


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
    '--fmt', required=False, default='json',
    help='Output format: genbank | fasta | json')
@click.option(
    '--source', required=False, default='ncbi',
    help='The data source.')
@click.option(
    '--ids', required=True,
    help='A list of IDs, e.g. GenBank accession IDs.')
@click.option(
    '--stdout', required=False, default=False, is_flag=True,
    help='Write to stdout instead of file.')
@click.argument('out', type=click.Path())
@click.command()
def load(ids, out, batch, email, db, fmt, source, stdout):
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
        eprint('Loading data from NCBI.')
        RETMAX = 10**9
        accessions = read_accessions(ids)
        dbase = db
        Entrez.email = email
        batchsize = batch

        eprint('Batch size:', batchsize)
        eprint('Starting download:')

        with click.progressbar(  # stackoverlow, 3173320
            accessions_to_fmt(
                accessions, dbase, batchsize, RETMAX, fmt),
                length=len(accessions),
                file=sys.stderr, show_pos=True, label='Progress:',
                fill_char=':', empty_char='.'
                ) as bar:
            for acc, record in bar:
                if not stdout:
                    write_record(out, acc, record, fmt)
                else:
                    print(record)  # This goes to stdout and can be piped.
        click.echo('Done.', err=True)


@click.option('--client', required=False, default='localhost:27017')
@click.option('--db', required=False)
@click.option('--cell', required=False, help='Cell name.')
@click.option('--client_env', envvar='ZOOCLIENT', required=False, default=None)
@click.option('--db_env', envvar='ZOODB', required=False, default=None)
@click.option('--cell_env', envvar='ZOOCELL', required=False, default=None)
@click.option('--delim', default='|')
@click.option('--fmt', required=True, default='json')
@click.option(
    '--query', required=False, default=None,
    help='A JSON file with the (pymongo) query syntax.')
@click.option(
    '--selection', default=None,
    help='Comma-seperated list of attributes to use in header.')
@click.argument('file', type=click.File('w+'))
@click.command()
def dump(
        client, client_env, db, db_env, cell, cell_env,
        query, selection, delim, file, fmt):
    '''Dump a query (maybe restricted by selection) to <fmt>.

    Example:

    \b
    cat q.json
    # {}
    zoo dump --query q.json --selection _id,meta.date,meta.geo.cou,seq \
    --delim "|" --fmt fasta dump.fa
    head -n2 dump.fa
    # >a0b5d956-a940-427d-b5ff-f3a22e750389|2015-09-07|sierra_leone
    # NNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTTAGGATCTTTTGTGTGCGAATAACTAT...

    \b
    zoo dump --query q.json --selection _id - > qdump.json
    zoo dump --query q.json --selection _id | wc -l
    zoo dump --query q.json --selection _id,seq --fmt fasta - | grep ">" | head
    # minhash
    zoo dump --query q.json --selection _id,seq --fmt fasta - | \
    sourmash compute -k 16 -n 100 --singleton --out q.sig -
    '''

    c = env_switch(client, client_env, db, db_env, cell, cell_env)
    records = eval_query(c, query, selection)  # will return selection as dict

    if fmt == 'json':
        for rec in records:
            try:
                file.write(json.dumps(rec) + '\n')
            except BrokenPipeError:  # somehow occurs when piping into head
                return

    elif fmt == 'fasta':
        s = selection.split(',')
        assert 'seq' in s  # a valid fasta needs a sequence
        for rec in records:
            header = delim.join(
                [deep_get(rec, key) for key in s if key != 'seq'])
            lines = '>{}\n{}\n'.format(header, rec['seq'])
            try:
                file.write(lines)
            except BrokenPipeError:
                return

