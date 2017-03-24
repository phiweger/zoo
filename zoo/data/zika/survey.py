from bson.json_util import dumps  # stackoverflow, 30333299
from copy import deepcopy
import json
from pyfaidx import Fasta
from pymongo import MongoClient
from uuid import uuid4
import zoo


client = MongoClient('localhost:27017')
db = client['zika']


with open(zoo.get_schema('base.json')) as infile:
    base = json.load(infile)


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

with open('PRJNA344504.json', 'w+') as outjson:
    outjson.write(dumps(db.survey.find()))
# fasta is 350 kb, JSON 359
