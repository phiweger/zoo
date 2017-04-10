import json
from jsonschema import validate
from pyfaidx import Fasta
from pymongo import MongoClient
from zoo import get_schema


def link_fasta(link, func=None):
    '''Given an identifier (e.g. UUID), access seq in supplement fasta file.

    func .. a key function to pass to pyfaidx.Fasta(), e.g.
    func = lambda x: x.split('|')[0]

    The access is random, i.e. quick, given that an index exists. If no index
    is found (file.fai) one is created.
    '''
    fa = Fasta(link.target, key_function=func)
    return str(fa[link.key])


def link_fastq():
    '''TODO'''
    pass


class Link():
    '''A link object.

    d = {'k': 'ebola.makona._id', 'v': 'a0b5d956-a940-427d-b5ff-f3a22e750389'}
    l = Link(d['k'], d['v'])
    if l.is_valid():
        print('Link is valid.')
    deep_get(next(l.access_document()), 'meta.geo.cou')
    'sierra_leone'

    l = Link(get_data('plum/NC_001445.fa'), 'NC_001445.1')
    l.access_fasta(func=lambda x: x.split(' ')[0])[:10]
    # 'AAAATATAAA'
    '''
    def __init__(self, target, key, **kwargs):
        '''
        supported target types: field ("internal" link), fasta
        '''
        self.target = target
        self.key = key

    def is_valid(self, schema=None):
        '''Validate link against a link schema.'''
        if not schema:
            with open(get_schema('fragments/link.json'), 'r+') as file:
                schema = json.load(file)
        assert validate(self.__dict__, schema=schema) is None
        return True

    def __repr__(self):
        gen = ('{}: {}\n'.format(k, v) for k, v in self.__dict__.items())
        out = ''
        for i in gen:
            out += i
        return out

    def access_document(self, client='localhost:27017'):
        u = self.target.split('.')
        db, cell = u[0], u[1]
        field = '.'.join(u[2:])
        c = MongoClient(client)[db][cell]
        return c.find({field: self.key})  # returns > 1 docs if key not unique

    def access_fasta(self, func=None):
        try:
            return link_fasta(self, func)
        except KeyError:
            print('The key indicated in the link is not found in file header.')
