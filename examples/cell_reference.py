# from itertools import islice
from bson.json_util import dumps
from copy import deepcopy
import json
import progressbar
from pyfaidx import Fasta
from pymongo import MongoClient
from uuid import uuid4
from zoo import get_schema
from zoo.utils import deep_get, deep_set


client = MongoClient('localhost:27017')
db = client['ref']


fp = 'viral.fna'
fa = Fasta(fp, key_function=lambda x: x.split('|')[1])


with open(get_schema('base.json')) as infile:
    schema = json.load(infile)


for i in fa:
    r = deepcopy(schema)
    r['_id'] = str(uuid4())
    r['sequence'] = str(i)
    u, v = i.name.split('.')
    deep_set(r, 'metadata.alt_id', {'gb': u, 'gb_version': v}, force=True)
    deep_set(
        r,
        'metadata.description',
        i.long_name.split('|')[-1].strip(),
        force=True)
    db.ref.insert_one(r)


with open('ref.json', 'w+') as outjson:
    outjson.write(dumps(db.ref.find(), indent=4))


from sourmash_lib import Estimators
from sourmash_lib.sbt import SBT, GraphFactory
from sourmash_lib.sbtmh import SigLeaf, search_minhashes
from sourmash_lib.signature import SourmashSignature

KSIZE = 16
N = 1000


# init SBT
factory = GraphFactory(ksize=KSIZE, starting_size=1e5, n_tables=4)
# 4 .. nt?
tree = SBT(factory, d=2)  # d .. see "n-ary " in notebook


bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
cursor = db.ref.find()
c = 0
for i in cursor:
    key = deep_get(i, 'metadata.alt_id.gb')
    seq = i['sequence']  # db.ref.find_one()['sequence']  # 'ACTG...'
    e = Estimators(ksize=KSIZE, n=N)
    e.add_sequence(seq, force=True)  # e.get_hashes()
    s = SourmashSignature(
        email='', estimator=e,
        name=key)

    leaf = SigLeaf(metadata=key, data=s)
    tree.add_node(node=leaf)
    c += 1
    bar.update(c)
# \ 9158 Elapsed Time: 0:01:49

# search the last fasta entry against the SBT (">0.95")
# filtered = tree.find(search_minhashes, s, 0.1)
# matches = [(str(i.metadata), i.data.similarity(s)) for i in filtered]
# [('0.95', 1.0)]  # fasta header, similarity


tree.save('ref')


'''
sourmash sbt_search -k 16 ref ~/repos/zoo/zoo/data/zika/survey.sig
# running sourmash subcommand: sbt_search
loaded query: survey... (k=16, DNA)
0.11 NC_012532
'''


record = next(db.ref.find({'metadata.alt_id.gb': 'NC_012532'}))
deep_get(record, 'metadata.description')
# 'Zika virus, complete genome'


# TODO: do the search internally, i.e. not via commanline




