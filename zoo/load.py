from Bio import Entrez, SeqIO
from io import StringIO
import json
import os
import sys
import urllib
from zoo.format import seqrecord2jsondict


def read_accessions(fp):
    with open(fp) as acc_lines:
        return [line.strip() for line in acc_lines]


def extract_records_gb(records_handle):
    '''...'''
    buffer = []
    for line in records_handle:
        if line.startswith('LOCUS') and buffer:
            # yield accession number and record
            yield buffer[0].split()[1], ''.join(buffer)
            buffer = [line]
        else:
            buffer.append(line)
    yield buffer[0].split()[1], ''.join(buffer)


def batch(accessions, size):
    '''Splitting accession list into batches of <batchsize>.'''
    l = len(accessions)
    for start in range(0, l, size):
        yield accessions[start:min(start + size, l)]


def process_batch(accessions_batch, db, batchsize, retmax, fmt):
    '''Process a batch of genbank accession IDs'''
    # get GI for query accessions
    query = ' '.join(accessions_batch)

    for attempt in range(10):  # stackoverflow, 2083987
        try:
            query_handle = Entrez.esearch(db=db, term=query, retmax=retmax)
            gi_list = Entrez.read(query_handle)['IdList']

            # get GB files
            search_handle = Entrez.epost(db=db, id=','.join(gi_list))
            try:
                search_results = Entrez.read(search_handle)
            except RuntimeError:
                print('\nMismatch btw/ selected database and accession IDs.')
                print('Aborted!')
                sys.exit()
            webenv = search_results['WebEnv']
            query_key = search_results['QueryKey']

            records_handle = Entrez.efetch(
                db=db, rettype='gb', retmax=batchsize,
                webenv=webenv, query_key=query_key)
        except (ConnectionResetError, TimeoutError, urllib.error.URLError) as e:
            print('\n')
            print(e)
            print('Reconnect, attempt', str(attempt + 1) + '/10.')
        else:
            break

    # divert here: fasta, gb or json?
    # Bio.Seq import record
    if fmt == 'genbank':
        yield from extract_records_gb(records_handle)
        # yield from, stackoverflow, 9708902
    elif fmt == 'json':
        for acc, r in extract_records_gb(records_handle):
            record = SeqIO.read(StringIO(r), 'genbank')
            yield acc, json.dumps(seqrecord2jsondict(record))
    elif fmt == 'fasta':
        for acc, r in extract_records_gb(records_handle):
            record = SeqIO.read(StringIO(r), 'genbank')
            yield acc, record.format('fasta')
    else:
        raise AttributeError(
                'Output format not supported, try fasta | genbank | json.'
                )


def accessions_to_fmt(accessions, db, batchsize, retmax, fmt):
    accession_batches = batch(accessions, batchsize)
    for acc_batch in accession_batches:
        yield from process_batch(acc_batch, db, batchsize, retmax, fmt)


def write_record(out, accession, record, fmt):
    ext = {'genbank': '.gb', 'fasta': '.fa'}
    if fmt in ['genbank', 'fasta']:
        if not os.path.exists(out):
            os.makedirs(out)
        with open(os.path.join(out, accession + ext[fmt]), 'w') as output:
            print(record, file=output)
    elif fmt == 'json':
        record = record.replace('sequence', 'seq')  # TODO: this is too hacky
        with open(out, 'a+') as output:
            print(record, file=output)
        # 'a+' append to file


