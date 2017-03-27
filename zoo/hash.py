from hashlib import md5
import json


def hash_document(d):
    '''md5 hash the JSON string representation of a dict w/o key "_id"'''
    _id = d.pop('_id')
    h = md5()
    h.update(json.dumps(d, sort_keys=True).encode('utf-8'))
    value = h.hexdigest()
    return _id, value
