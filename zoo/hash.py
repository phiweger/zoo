from hashlib import md5
from zoo.utils import ordered, flat


def hash_document(d):
    '''md5 hash the JSON string representation of a dict w/o key "_id"

    We need to be careful about the order of the nested dicts (i.e. JSON)
    entries when hashing. A different order would produce a different hash
    although the content is not changed. However, the two key objects
    Python uses when representing a JSON file, i.e. dicts and lists, do not
    enforce any order. And usually we do not care, so implementing everything
    with heavy objects such as OrderedDict seems overkill.

    When zoo computes a document's hash, it firsts recursively orders
    all the keys and values of the document and then flattens it. Note that
    ordered does sort tuples as well. They are not permitted in Python syntax,
    but it seems to make thinks more robust that way.
    '''
    _id = d.pop('_id')
    try:
        _ = d.pop('md5')  # discard old hash
    except KeyError:
        pass
    h = md5()
    f = flat(ordered(d))
    for i in f:
        h.update(str(i).encode('utf-8'))
    value = h.hexdigest()
    return _id, value



