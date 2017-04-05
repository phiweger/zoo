import hashlib
import json
import numpy as np
from pyfaidx import Fasta
from zoo import get_schema


def msa_subset(msa, X, y):
    '''
    subset MSA but DO KEEP X - y mapping

    Usage:
    a = msa_subset(msa, X_train, y_train)
    '''
    mp = {}

    for key, value in zip(X, y):
        mp[key] = [value]

    for i in msa:
        try:
            mp[i.metadata['id']].append(str(i))
        except KeyError:
            continue
    # example dict entry: 'KP285028': ['Avian', 'ATGGAGGACTTTGTGCGACAATGCT...]

    ids = []
    tag = []
    seq = []
    for key, values in mp.items():
        ids.append(key)
        tag.append(values[0])
        seq.append(list(values[1]))

    # return np.array(seq), np.array(ids), np.array(tag)
    return np.array(seq), np.array(ids), np.array(tag)


def align_deopt(reference, query):
    '''
    Given a sequence from a Mafft MSA as reference (usually wilde type),
    "align" a modified version of the reference (i.e. query), by putting
    in the gaps introduced into the reference by the MSA.

    We need this function to evaluate modified sequences with information
    obtained from an MSA based analysis (e.g. entropy, feature importance).

    Example:

    - reference:        ACTG
    - query:            AATG
    - reference in MSA: A-C--TG

    align_deopt(reference, query)
    # A-A--TG
    '''
    s = ''
    counter = 0
    for i in range(len(reference)):
        if reference[i] == '-':
            s += '-'
        else:
            s += query[counter]
            counter += 1
    assert s.replace('-', '') == query  # True
    return s


def encode_gaps(seq):
    '''
    Encode the gaps of a given sequence from an MSA e.g. to store together
    with an MSA identifier as a sequence derivative, for later MSA
    reconstruction.

    encode_gaps('A--CTGA---GGTAGGT-AA')
    # {1: 2, 7: 3, 17: 1}
    # more efficient than list of tuples like [(1, 2), ...]
    # stackoverflow, 15641344

    # alternative:
    # 10011110001111111011 compressed (hard to share bc/ binary?)
    # alternative:
    # list of all positions where gap
    '''
    bv = [1 if nt == '-' else 0 for nt in seq]  # b .. binary vector

    # count all progressions of 1, write to dict when terminated by 0
    d = {}
    index = 0
    counter = 0
    for i in bv:
        if i == 1:
            counter += 1
        else:
            if counter != 0:
                d[index - counter] = counter
                counter = 0  # reset to 0
        index += 1
    return d


def decode_gaps(seq, gapdict, uppercase=False):
    '''
    Decode the gaps of a given sequence to make it "fit" into an MSA.
    The default behavior is to not meddle with the upper/ lower case.

    decode_gaps('ACTGAGGTAGGTAA', {1: 2, 7: 3, 17: 1})
    # 'A--CTGA---GGTAGGT-AA'
    '''
    result = ''
    index = 0
    for nt in seq:
        try:
            v = gapdict[index]
            insert = v * '-'
            result = result + insert + nt
            index += len(nt + insert)
        except KeyError:
            result += nt
            index += 1
    return result


def hash_seq(seq_iterator, method='md5'):
    '''
    Given any sequence iterator, calculate its hash with <method>. We use this
    to hash an MSA to an ID as well as verify its integrity (think checksum).

    Note that order does matter for most hashing methods, so internally the
    sequence iterator is lexicographically sorted before fed to the hash func.

    Examples:

    hash_seq(['AAAA', 'ACTG'])
    # '3b0d5a673fb29e5201e4587a35e2576d'

    hash_seq(['ACTG', 'AAAA']) == hash_seq(['AAAA', 'ACTG'])
    # True

    from pyfaidx import Fasta
    fp = 'msa.mafft.fa'
    fa = (str(i).upper() for i in Fasta(fp))  # Mafft output is lowercase
    hash_seq(fa)
    # 'bba3a154dae681c96...'
    '''
    m = getattr(hashlib, method)()  # stackoverflow, 3061
    for i in sorted(seq_iterator):  # order does matter for md5
        m.update(i.encode('utf-8'))
    return m.hexdigest()


# think of msas accross multiple collections
def import_msa(db, fp, label=''):
    ''' Import MSA members into respective docs, use gap notation.

    We want to preserve a given MSA while not storing it as a flat file
    on disk (an forgetting after 2 days which sequences were in it etc.)
    Thus, the gaps of each member sequence of the (unconstrained) MSA
    are encoded end saved to the respective document under
    "derivative.msa". The full MSA can be reconstructed from the MSA id.
    Note that we can reconstruct a part of the MSA just as well by restricting
    our MSA ID query to only a subset of the collections whose documents
    went into the MSA.

    We assume that the MSA has the form:

    >UUID1|collectionA
    act--atga--AAGcc...
    >UUID2|collectionA
    more--se-quence-...
    >UUID3|collectionB
    even--more------...

    This can be achieved with zoo.utils.to_fasta().

    Note that the case of the nucleotide seq is not touched by zoo. If it
    is stored in uppercase in zoo, it is exported that way to a format
    as input to some MSA algo. Hashing the MSA is case sensitive. Usually
    this should not matter, but if you experience any MSA ID discrepancy,
    this is a likely source.


    client = MongoClient("localhost:27017")
    db = client["jmtv"]
    '''
    with open(get_schema('derivative_msa.json')) as infile:
        schema = json.load(infile)

    msa = Fasta(fp)
    gen = (str(i).upper() for i in msa)
    h = hash_seq(gen)

    for i in msa:
        name, collection = i.name.split('|')
        entry = schema.copy()  # copy() necessary here?
        entry['id'] = h
        entry['label'] = label
        entry['gaps'] = encode_gaps(str(i))

        c = db.get_collection(collection)
        c.update_one(
            {'_id': name},
            {'$push': {'derivative.msa': entry}}
            )
        # https://docs.mongodb.com/manual/reference/operator/update/push/


# def export_msa():
#   pass

















