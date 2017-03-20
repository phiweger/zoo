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


'''
get a number of docs, to SBT in tempfile or disk, query an interesting
seq, save result in doc?

sbt_index([cursor], fp='', temporary=True)
sbt_combine([list, of, cursors], fp='', temporary=True)
sbt_search(sbt, query)
'''


# https://github.com/luizirber/2016-sbt-minhash/blob/master/notebooks/SBT%20with%20MinHash%20leaves.ipynb
# http://blog.luizirber.org/2016/12/28/soursigs-arch-1/
from sourmash_lib import Estimators
from sourmash_lib.sbt import SBT, GraphFactory
from sourmash_lib.sbtmh import SigLeaf, search_minhashes
from sourmash_lib.signature import SourmashSignature

seq = db.ref.find_one()['sequence']  # 'ACTG...'
e = Estimators(ksize=16, n=200)
e.add_sequence(seq)  # e.get_hashes()
s = SourmashSignature('foo', e, name='')
factory = GraphFactory(ksize=16, starting_size=1e5, n_tables=4)
# 4 .. nt?
tree = SBT(factory)  # d .. see "n-ary " in notebook
'''
How does the internal nodes / total ratio affect query times? Test and
put as section in master thesis.
'''

leaf = SigLeaf(metadata='bar', data=s)
# SigLeaf(metadata, data, name=None)

tree.add_node(node=leaf)
tree.print()  # ignore pylint
filtered = tree.find(search_minhashes, s, 0.1)

matches = [(str(i.metadata), i.data.similarity(s)) for i in filtered]
# seems to work
# [('bar', 1.0)]










