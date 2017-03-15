import numpy as np

# DEPRECATED
# def msa_subset(msa, ids):
#     '''
#     Returns the sequences whose id is in ids.

#     Usage:
#     a = msa_subset(msa, list(X_train))
#     next(a)
#     # 'ATGAAT...'

#     '''
#     for line in msa:
#         if line.metadata['id'] in ids:
#             yield str(line)


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
    index = 0
    result = ''
    for nt in seq:
        try:
            n = gapdict[index]
            result += n * '-' + nt
            index = index + n
        except KeyError:
            result += nt
        index += 1
    if uppercase is True:
        result = result.upper()
    return result
