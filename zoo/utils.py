'''
Careful not to mix up methods between zoo and sourmash.

TODO: iter_seq: customize name like _id|host|...
'''


from __future__ import print_function  # Python 2 compatible
from functools import reduce
import hashlib
import json
import os
import progressbar
from pyfaidx import Fasta
import re
import sourmash_lib
from sourmash_lib import signature
import sys
import tempfile


def eprint(*args, **kwargs):
    '''Print to stderr, stackoverflow, 5574702.'''
    print(*args, file=sys.stderr, **kwargs)


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

    Naming: name is constructed from header argument (a list), will
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


def get_id(collection, ids):
    '''
    In: a list of IDs and a mongoDB collection
    Out: a cursor of documents corresponding to the IDs

    q = get_id(collection, X_test)
    q.count()  # 3105
    len(X_test)  # 3105

    # just for illustration, as this will consume the cursor
    from collections import Counter
    Counter(i['metadata']['host'] for i in q_test)
    Counter({'Avian': 1035, 'Human': 1035, 'Swine': 1035})
    '''
    match = {'_id': {'$in': list(ids)}}
    q = collection.find(match)
    return q


def deep_get(obj, item):  # , fallback=None  # see comment below
    '''Steps through an item chain to get the ultimate value.

    Implementation inspired by: stackoverflow, 25833613, alternative:
    www.haykranen.nl/2016/02/13/handling-complex-nested-dicts-in-python/

    If ultimate value or path to value does not exist, does not raise
    an exception and instead returns `fallback`.

    d = {'snl_final': {'about': {'_icsd': {'icsd_id': 1}}}}
    deep_get(d, 'snl_final.about._icsd.icsd_id')
    # 1
    deep_get(d, 'snl_final.about._sandbox.sbx_id')
    # (None returned)
    '''

    # Use this when no error raising is desired. Can't distinguish
    # btw/ KeyError and {'key_present': None}. This is ennoying when
    # working with nested dicts, so commented out, for reference only.
    # def getitem(obj, name):
    #     try:
    #         return obj[name]
    #     except (KeyError, TypeError):
    #         # KeyError: 'k' not present
    #         # TypeError: 'k.a' not present, "not subscriptable"
    #         return fallback

    def getitem(obj, name):
        return obj[name]  # err out if key does not exist

    return reduce(getitem, item.split('.'), obj)


def deep_set(d, key, value, force=False, replace=False):
    '''Dynamically sets values deep in a dict.

    If keys do not exist, they are created. Note the divergence from the
    original implementation, bc/ did not see the point for line "dd = d" in
    stackoverflow, 21297475.

    Example:

    d = {}
    deep_set(d, 'a.b.c', 'X', force=True)
    deep_set(d, 'a.b.d', 'Y', force=True)
    deep_set(d, 'a.b.d', 'Y2')  # not allowed
    deep_set(d, 'a.b.d', 'Y2', replace=True) # allowed
    deep_set(d, 'a.b.d', {'foo': 'bar'})

    # append to list, pop from dict etc.
    deep_get(d, 'a.b.d').pop('foo')
    # 'bar'
    # Note: "get" not "set" is used here.

    This function took quite some time to figure out, so left here for
    future reference:

    - in
        https://github.com/mahmoud/boltons
      see
        http://sedimental.org/remap.html
    - looks promising: stackoverflow, 3797957
    - stackoverflow, 3031219: from dotDict import dotdictify
    - stackoverflow, 3232943: Alex Martelli answer
    - stackoverflow, 11918852
    '''
    keys = key.split('.')
    latest = keys.pop()

    # see deep_get comments on fallback changes
    # key_exists = deep_get(d, key)  # , fallback=None
    # if key_exists:

    try:
        _ = deep_get(d, key)
        key_exists = True
    except (KeyError, TypeError):
        key_exists = False
    # deep_get erring out means key not there, no error - key present,
    # this can distinguish KeyError and {'key_present': None} in dict
    if key_exists:
        if replace is False and _:
            # i.e. field not None, "", [], ...
            # i.e. field is non-False and we don't want to replace it
            raise TypeError(
                'Key exists, item assignment not allowed w/ replace=False')
        else:
            reduce(dict.__getitem__, keys, d)[latest] = value
            # stackoverflow, 11918852

    # Else, create a new key if force=True. Otherwise err out.
    else:
        if not force:
            raise KeyError(
                'Key not present. Use "force=True" to create key.')
        else:
            for k in keys:
                d = d.setdefault(k, {})
            d.setdefault(latest, value)


def ordered(obj):
    '''Sort a nested dict's keys AND values recursively.

    stackoverflow, 25851183

    Example:

    import json

    a = json.loads("""
    {
        "errors": [
            {"error": "invalid", "field": "email"},
            {"field": "name", "error": "required", "lala": {"new": 4}}
        ],
        "list": [5,3,2,4,1,6],
        "success": false
    }
    """)

    b = json.loads("""
    {
        "success": false,
        "errors": [
            {"lala": {"new": 4}, "error": "required", "field": "name"},
            {"error": "invalid", "field": "email"}
        ],
        "list": [1,2,3,4,5,6]
    }
    """)

    assert ordered(a) == ordered(b)
    '''
    if isinstance(obj, dict):
        return sorted((k, ordered(v)) for k, v in obj.items())
    if isinstance(obj, list):
        return sorted(ordered(x) for x in obj)
    else:
        return obj


def flat(l):
    '''
    Given a nested list and/ or tuple flatten (= "unnest") it.
    '''
    result = []
    for item in l:
        if isinstance(item, (list, tuple)):
            # http://bit.ly/2biB44i
            # http://bit.ly/2b3Nwa4
            result.extend(flat(item))
        else:
            result.append(item)
    return result
