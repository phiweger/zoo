# from itertools import islice
import json
import progressbar
from pyfaidx import Fasta
from pymongo import MongoClient
from uuid import uuid4
from zoo.dotmap import DotMap
from zoo import get_schema


client = MongoClient('localhost:27017')
db = client['mock']


fp = '/Users/pi/data/reference/virus/viral.1.1.genomic.fna'
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
