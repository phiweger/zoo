'''
Careful not to mix up methods between zoo and sourmash.

TODO: iter_seq: customize name like _id|host|...
'''

import hashlib
import json
import re
import os
import progressbar
import sourmash_lib
from sourmash_lib import signature
import tempfile


def parse_iav(isolate):
    '''
    parse the data from influenza name into dict

    test:
    a = 'A/mallard/Interior Alaska/6MP0155/2006(H3N8)'
    b = 'A/blue-winged Teal/Minnesota/AI09-2977/2009(H4N8)'
    c = 'A/Peru/PER175/2012(H3N2)'
    '''

    # regexp structure check

    feat = isolate.split('/')
    assert feat[0] == 'A', 'Not Influenza A Virus.'

    subtype = re.search(r'\((.*?)\)', isolate).group(1)
    feat[-1] = feat[-1].split('(')[0]
    if len(feat) == 4:
        host = 'human'
        geo = feat[1]
    else:
        host = feat[1].lower()
        geo = feat[2]

    return({
        'influenza': feat[0],
        'host': host,
        'subtype': subtype,
        'geolocation': geo
        })


def parse_date(date):
    '''
    Parse (influenza) dates (source: NCBI Genomes, i.e. GenBank) into
    nicely searchable format for insertion in zoo database. Four formats
    are recognized:

    - 1976
    - 2011/01
    - 2011/01/27  # ymd
    - ('NON', 'Unknown', '')

    Usage:

    parse_date('2019/08')
    # {'d': '', 'm': 8, 'y': 2019}
    '''
    record = {'y': '', 'm': '', 'd': ''}
    if date in ('NON', 'Unknown', 'unknown', ''):
        return record
    else:
        s = date.split('/')
        for i in zip(range(3), 'y m d'.split(' ')):
            j, tag = i
            try:
                record[tag] = int(s[j])
            except IndexError:
                continue
        return record


def locate(locations, seq):
    '''
    Return the (concatenated) sequence string corresponding to the
    list of locations (a location being a tuple of start and end
    coordinates).

    sample.annotation.M1.location
    # ['6-764']

    764 - 5
    # 759

    # len(locate(sample.annotation.M1.location, sample.sequence))
    759

    # checked that it works for 2 locations as well, should be do so for > 2
    '''
    import re

    product = ''
    for location in locations:
        start, end = re.findall('\\d+', location)
        fragment = seq[int(start) - 1:int(end)]
        # -1 because GenBank entries start indexing at 1, not 0 as in python
        # when slicing e.g. [1:5], python excludes position 5, so here no -1
        product += fragment
        # order matters, eyeballed some entries and they seem ordered
    return product


def to_fasta(fp, iterable):
    '''
    Write what comes out of iter_seq to fasta
    '''
    with open(fp, 'w+') as outfile:
        for i in iterable:
            header, seq = i
            outfile.write(
                '>{}\n{}\n'.format(
                    header, seq
                    )
                )


def iter_seq(query, annotation=None, field_seq='sequence', header=[]):
    '''
    > id
    seq

    hosts = 'Human Avian Swine'.split(' ')
    query = col.find({
        'annotation.name': 'HA',
        'metadata.host': {'$in': hosts},
        'metadata.id': {'$nin': ['22271']}
    })

    a = iter_seq(annotation='HA', query=query)
    next(a)
    # ('CY196419',
    # 'ATGATGGTCA...

    Note that a tuple is returned: (_id, sequence)

    If no annotation is specified, the raw sequence is returned.

    naming: name is constructed from header argument (a list), will
    concatenate with "|", for this info the 'metadata' field is queried,
    and the '_id' field is appended as well, both hardcoded, deal
    '''
    if annotation is None:
        for doc in query:
            name = doc['_id']
            seq = doc[field_seq]
            yield name, seq
    else:
        for doc in query:
            name = doc['_id']
            if header:
                name += '|' + '|'.join([doc['metadata'][h] for h in header])

            for i in doc['annotation']:
                if i['name'] == annotation:
                    seq = locate(i['location'], doc[field_seq])
                    yield name, seq  # or customize _id|host|...


def minhash(sequence_gen, ksize, n, **kwargs):
    '''
    as kwargs we can pass (=default)
    - with_cardinality=False
    - track_abundance=False

    Note that a tuple is returned: (_id, estimator), i.e. we continue
    to carry the _id along.
    '''
    for name, seq in sequence_gen:
        e = sourmash_lib.Estimators(
            n=n,
            ksize=ksize,
            **kwargs
            )
        e.add_sequence(seq, force=True)
        yield name, e


def save_minhash(minhash_gen, handle=None, email=''):  #
    '''
    issue:
    https://github.com/dib-lab/sourmash/issues/131
    suggested lead:
    https://github.com/dib-lab/sourmash/blob/master/utils/compute-dna-mh-another-way.py
    relevant set of functions:
    https://github.com/dib-lab/sourmash/blob/master/sourmash_lib/signature.py

    from itertools import islice
    print(save_minhash(islice(gen, 2)))

    fp = '/some/path/to.json'
    with open(fp, 'w+') as outfile:
        save_minhash(islice(gen, 2), handle=outfile)
    '''
    l = []
    bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
    counter = 0

    print('Generating signatures ...')
    for mh in minhash_gen:
        name, e = mh
        s = signature.SourmashSignature(email, e, name=name, filename='zoo')
        l.append(s)
        # load all to memory ... bad, TODO: stream into file handle

        counter += 1
        bar.update(counter)

    print('\nSaving signatures ...')
    if handle is None:
        return signature.save_signatures(l)  # [s] instead of l
    else:
        return signature.save_signatures(l, handle)


def create_sbt(seq_iterable):
    '''
    Given a query against mongodb, create a SBT from them.

    Allow SBT construction of arbitrary sequences, like RdRp

    input any sequence iterable (ideally generator)
    '''
    pass


def parse_categorize(csv, log):
    '''
    parse_categorize(fp + name_query + '.csv', fp + name_query + '.log')
    '''
    with open(csv, 'r+') as incsv, open(log, 'r+') as inlog:
        for line in inlog:
            if 'for' in line:
                a = line.split(' ')[1].strip(',')
                b, c = incsv.readline().strip().split(',')[1:]
                yield a, b, c


def minhash_categorize(fp, name_query, name_sbt):
    '''
    convenience function to call sourmash categorize from within
    python and parse the awkwardly formatted result automatically into
    a nice result.out

    in: pointers to two files, all calculations performed outside of
    python session

    out: result.out

    equivalent to:
    sourmash categorize -k 16 --csv control.csv ... 2> file.log
    '''
    returncode = os.system(
        ' '.join([
            'sourmash categorize',
            '--csv', fp + name_query + '.csv',
            fp + name_sbt,
            fp + name_query + '.json',
            '2>&1 | tee ' + fp + name_query + '.log'
            ])
        # redirect stderr to stdout (the pipe) and into tee, which
        # prints it and saves it at the same time
        )

    if returncode == 0:
        print(
            'sourmash categorize successful, parsing result and saving to:\n',
            fp + name_query + '.out'
            )

    with open(fp + name_query + '.categorize', 'w+') as outfile:
        a = parse_categorize(
            fp + name_query + '.csv', fp + name_query + '.log'
            )
        for i in a:
            outfile.write('\t'.join(i))
            outfile.write('\n')


def minhash_search(fp, name_query, name_sbt, k, threshold=0.08):
    '''
    usage:
    minhash_search(fp, name_query, name_sbt, threshold=0.7)


    Note one important thing: the "name" field in the query.json
    cannot contain any characters that are not
    '''

    # load json
    with open(fp + name_query + '.json', 'r+') as infile:
        file = json.load(infile)
    # create new dir
    os.system('mkdir ' + fp + name_query)

    # save things to tempdir
    with tempfile.TemporaryDirectory() as tempdir:
        for i in file:
            # with tempfile.NamedTemporaryFile('w+') as temp:
            filename = ''.join([tempdir, '/', i['name'], '.json'])
            # we now hash filename to avoid OS file name problems, e.g. "|"
            h = hashlib.md5()
            h.update(filename.encode('utf-8'))
            filehash = h.hexdigest()
            fp_out = fp + name_query + '/' + filehash

            with open(fp_out, 'w+') as outfile:
                outfile.write('# ' + i['name'] + '\n')

            with open(filehash, 'w+') as temp:
                temp.write('[')  # sourmash sbt_search requires these brackets
                json.dump(i, temp, indent=True)
                temp.write(']')
                temp.seek(0)  # this line cost me 1 h
                # print(temp.name)
                # temp is closed as soon as we leave the (with ...) context
                command = ' '.join([
                    'sourmash sbt_search',
                    '-k', str(k),
                    '--threshold', str(threshold),
                    fp + name_sbt + '.sbt.json',
                    filehash,
                    # '> ' + fp + name_query + '/' + i['name'] + '.search'
                    '>> ' + fp_out
                    ])
                os.system(command)

# def select_taxonomy(id):
#     '''
#     Given an id, return all entries from all collections associated
#     with taxonomy id.

#     Aim at some point is to build an RNA-virus SBT and query it.
#     '''
#     pass




# def export_fasta(cursor):
#     '''
#     Given a mongodb cursor (i.e. a query), export the sequence to fasta
#     w/ the indicated fields concatenated by "|".
#     '''
#     pass


# # def minhash(cursor):
# #     '''
# #     Given a document set (a query), calculate minhash sequence with
# #     k = {15, 30} (sourmash).
# #     '''
# #     pass


# def sbt_index(cursor):
#     '''
#     Given a query, use sourmash to construct an SBT .. sequence bloom tree.
#     '''
#     pass


# def phylo_query(signature, sbt, n):
#     '''
#     Given a minhash signature calculated by sourmash, return the n "closest"
#     documents.
#     '''
#     pass


# def complete_iav(cursor):
#     '''
#     Return a generator w/ entries that have complete genomes, i.e. entries for
#     8 segments.
#     '''
#     pass