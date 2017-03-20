# from itertools import islice
import json
import progressbar
from pyfaidx import Fasta
from pymongo import MongoClient
from uuid import uuid4
from zoo import get_schema
from zoo.dotmap import DotMap
from zoo.utils import deep_get


client = MongoClient('localhost:27017')
db = client['mock']


fp = '/path/to/viral.1.1.genomic.fna'
fa = Fasta(fp, key_function=lambda x: x.split('|')[1])


with open(get_schema('base.json')) as infile:
    schema = json.load(infile)


bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
count = 0
for i in fa:
    d = DotMap(schema)

    d._id = str(uuid4())
    d.sequence = str(i)
    d.metadata.ids = [i.name]
    # d.to_dict()
    # d.dget_attr('metadata.ncbi')
    db.ref.insert_one(d.to_dict())
    count += 1
    bar.update(count)

DotMap(db.ref.find_one()).dget_attr('metadata.ids')
# ['NC_021865.1']


def to_fasta(dmap, header_attr=[]):
    l = []
    for i in header_attr:
        l.append(dmap.dget_attr(i))
    print(l)
    print('|'.join(l))
    print(dmap.sequence[:10])


d = {'a': {'b': [1, 2, 3], 'c': {
    'd': 4}, 'd': [{'hello': 'A'}, {'hello': 'B'}]}, 'b': {}}

deep_get(d, 'a.d', fallback=None)
# [{'hello': 'A'}, {'hello': 'B'}]
deep_get(d, 'a.c.d.e.f.g', fallback='foo')
# 'foo'








