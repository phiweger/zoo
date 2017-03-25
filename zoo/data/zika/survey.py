from bson.json_util import dumps  # stackoverflow, 30333299
from copy import deepcopy
import json
from pyfaidx import Fasta
from pymongo import MongoClient
from sourmash_lib import Estimators, signature
from uuid import uuid4
import zoo


client = MongoClient('localhost:27017')
db = client['zika']


with open(zoo.get_schema('base.json')) as infile:
    base = json.load(infile)

fa = Fasta(
    zoo.get_data('zika/survey.fa'),
    key_function=lambda x: x.split('|')[1])

# load fasta
# >ENA|KY014295|KY014295.2 Zika virus isolate Zika virus/H.s...
fa = Fasta('survey.fa', key_function=lambda x: x.split('|')[1])

# write keys to txt
with open('survey.txt', 'w+') as outfile:
    for k in fa.keys():
        outfile.write(k + '\n')

# create fasta
with open('survey.fa', 'w+') as outfile:
    for i in fa:
        _id = str(uuid4())
        seq = str(i)

        # create db entry
        record = deepcopy(base)
        record['sequence'] = str(i)
        record['_id'] = _id
        db.survey.insert_one(record)

# db.survey.count()
# 33

with open('survey.json', 'w+') as outjson:
    outjson.write(dumps(db.survey.find(), indent=4))
# fasta is 350 kb, JSON 359


# now create a minhash of the entire sequences to search
e16 = Estimators(ksize=16, n=1000)
e31 = Estimators(ksize=31, n=1000)
fn = 'survey.sig'

cursor = db.survey.find({}, {'sequence': 1, '_id': 0})
for record in cursor:
    for k, v in record.items():
        e16.add_sequence(v, force=True)
        e31.add_sequence(v, force=True)
        # force bc/ ValueError: invalid DNA character in sequence: Y
s16 = signature.SourmashSignature(
    email='',
    estimator=e16,
    name='survey',
    filename=fn
    )
s31 = signature.SourmashSignature(
    email='',
    estimator=e31,
    name='survey',
    filename=fn
    )

with open(fn, 'w+') as outsig:
    signature.save_signatures([s16, s31], fp=outsig)
