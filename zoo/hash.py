import hashlib
from zoo.utils import ordered, flat


def hash_content(d, hashfun='sha256'):
    '''Hash the content of a hashmap without converns for order.

    We need to be careful about the order of the nested dicts (i.e. JSON)
    entries when hashing. A different order would produce a different hash
    although the content has not changed.
    When zoo computes a document's hash, it firsts recursively orders
    all the keys and values of the document and then flattens it. Note that
    ordered does sort tuples as well. They are not permitted in JSON syntax,
    but it seems to make thinks more robust that way.
    '''
    fun = getattr(hashlib, hashfun)
    h = fun()
    f = flat(ordered(d))
    for i in f:
        print(i)
        h.update(str(i).encode('utf-8'))
    return h.hexdigest()
