from pymongo import MongoClient
from uuid import uuid4
from zoo import get_data
from zoo.link import Link
from zoo.utils import deep_get


def test_link_document():
    c = MongoClient('localhost:27017')['test']['cell']
    _id = str(uuid4())
    c.insert_one({
        '_id': _id,
        'nested_key': {'key': 'brown'},
        'nested_value': {'value': 'fox'}
        })
    d = {'k': 'test.cell._id', 'v': _id}  # recover _id
    l = Link(d['k'], d['v'])
    assert l.is_valid()
    assert deep_get(next(l.access_document()), 'nested_value.value') == 'fox'


def test_link_document_nested_key():
    c = MongoClient('localhost:27017')['test']['cell']
    _id = str(uuid4())
    nested_key = str(uuid4())
    c.insert_one({
        '_id': _id,
        'nested_key': {'key': nested_key},
        'nested_value': {'value': 'fox'}
        })
    d = {'k': 'test.cell.nested_key.key', 'v': nested_key}
    l = Link(d['k'], d['v'])
    assert l.is_valid()
    assert deep_get(next(l.access_document()), 'nested_value.value') == 'fox'


def test_link_document_multiple_key_matches():
    c = MongoClient('localhost:27017')['test']['cell']
    nested_key = str(uuid4())
    c.insert_one({
        '_id': str(uuid4()),  # ID1
        'nested_key': {'key': nested_key},  # same nested key
        'nested_value': {'value': 'fox'}
        })
    c.insert_one({
        '_id': str(uuid4()),  # ID2
        'nested_key': {'key': nested_key},
        'nested_value': {'value': 'jumps'}
        })
    d = {'k': 'test.cell.nested_key.key', 'v': nested_key}
    l = Link(d['k'], d['v'])
    s = set()
    for i in l.access_document():
        s.update([deep_get(i, 'nested_value.value')])
    assert set(['fox', 'jumps']) == s


def test_link_fasta():
    l = Link(get_data('plum/plum.fa'), 'NC_014697.1')
    seq = l.access_fasta(func=lambda x: x.split(' ')[0])[:10]
    assert seq == 'TGGGCGAACG'
