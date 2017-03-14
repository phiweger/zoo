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
