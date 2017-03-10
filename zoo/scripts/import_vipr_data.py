import pymongo
from pymongo import MongoClient
import json
import operator
import progressbar
from pyfaidx import Fasta
import sys


sys.path.append('/Users/pi/repos/zoo')
from zoo.utils import parse_date


client = MongoClient("localhost:27017")
db = client["mock"]  # not zoo straight away
collection = db.get_collection('flavivirus')


client = MongoClient("localhost:27017")
db = client["zoo"]  # not zoo straight away
collection = db.get_collection('influenza_a_virus')


'''
{'_id': 'CY052137',
 'annotation': [{'accession': 'ACZ98217',
   'fuzzy': 0,
   'location': ['2-2281'],
   'name': 'PB2'}],
 'metadata': {'age': '14Y',
  'country': 'USA',
  'date': {'d': 8, 'm': 6, 'y': 2009},
  'gender': '',
  'host': 'Human',
  'id': '1005127',
  'isolate': 'Influenza A virus (A/New York/3976/2009(H1N1))',
  'segment_number': '1',
  'seqlength': '2293',
  'subtype': 'H1N1'},
 'sequence': 'TATGGAGAGAATAAAAG
'''


fp = '/Users/pi/data/vipr/test/Results.tsv'
with open(fp, 'r+') as infile:
    header = infile.readline()

header = header.strip().split('\t')
'''
['Strain Name',
 'Virus Type',
 'Subtype/Genotype (ViPR)',
 'GenBank Accession',
 'Sequence Length',
 'Collection Date',
 'Host',
 'GenBank Host',
 'Country',
 'Mol Type']
'''


'''
Proposal for virus document.
'''


# initialize document structure
d = dict.fromkeys('_id metadata taxonomy annotation sequence'.split(' '))
d['annotation'] = []
# an element in the annotation subdocument
anno = dict.fromkeys('accession fuzzy location name'.split(' '))
d['annotation'] = []
# the structure of the metadata
meta = dict.fromkeys(
    'sample misc sequence'.split(' '))
# misc contains virus-specific info, like segment number, subtype, ...
misc = dict.fromkeys('moltype host_gb'.split(' '))
# sample (= isolate) contains info about sample taken, minus time and location
sample = dict.fromkeys('country host date age id gender geo'.split(' '))

tax = dict.fromkeys('order genus family species strain other'.split(' '))
tax['other'] = []  # type, subtype, clade, cluster
seq = dict.fromkeys('length gc'.split(' '))
meta['misc'] = misc
meta['sample'] = sample
d['taxonomy'] = tax
meta['sequence'] = seq
d['metadata'] = meta
d['metadata']['sample']['geo'] = {'lat': None, 'long': None}
d['metadata']['sample']['date'] = {'y': None, 'm': None, 'd': None}


print(json.dumps(d, sort_keys=True, indent=4))
'''
{
    "_id": null,
    "annotation": [],
    "metadata": {
        "misc": {
            "host_gb": null,
            "moltype": null
        },
        "name": null,
        "sample": {
            "age": null,
            "country": null,
            "date": null,
            "gender": null,
            "host": null,
            "id": null
        },
        "sequence": {
            "gc": null,
            "length": null
        }
    },
    "sequence": null,
    "taxonomy": {
        "subtype": null,
        "type": null
    }
}
'''


fp = '/Users/pi/data/vipr/flaviviruses/metadata.tsv'
bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
counter = 0

with open(fp, 'r+') as infile:
    text = infile.readlines()
    text = text[1:]
    for line in text:
        counter += 1
        e = line.strip().split('\t')

        # fill document
        doc = d.copy()
        doc['taxonomy']['strain'] = e[0]
        doc['taxonomy']['type'] = e[1]
        doc['taxonomy']['subtype'] = e[2]
        doc['_id'] = e[3]
        doc['metadata']['sequence']['length'] = int(e[4])
        doc['metadata']['sample']['date'] = parse_date(e[5])
        doc['metadata']['sample']['host'] = e[6]
        doc['metadata']['misc']['host_gb'] = e[7]
        doc['metadata']['sample']['country'] = e[8]
        doc['metadata']['misc']['moltype'] = e[9]

        # print(json.dumps(doc, sort_keys=True, indent=4))
        try:
            db.flavivirus.insert_one(doc)
        except pymongo.errors.DuplicateKeyError:
            continue
        bar.update(counter)
'''
['Strain Name',
 'Virus Type',
 'Subtype/Genotype (ViPR)',
 'GenBank Accession',
 'Sequence Length',
 'Collection Date',
 'Host',
 'GenBank Host',
 'Country',
 'Mol Type']
'''

pipeline = [{'$group': {'_id': '$metadata.sample.host', 'cnt': {'$sum': 1}}}]
q = collection.aggregate(pipeline)

e = {}
for i in q:
    e[i['_id']] = i['cnt']
sorted_e = sorted(e.items(), key=operator.itemgetter(1))


'''
Coronavirus.
'''


fp = '/Users/pi/data/vipr/coronaviruses/metadata.tsv'
bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
counter = 0

with open(fp, 'r+') as infile:
    text = infile.readlines()
    text = text[1:]
    for line in text:
        counter += 1
        e = line.strip().split('\t')

        # fill document
        doc = d.copy()
        doc['taxonomy']['strain'] = e[0]
        doc['taxonomy']['species'] = e[1]
        doc['_id'] = e[2]
        doc['metadata']['sequence']['length'] = int(e[3])
        doc['metadata']['sample']['date'] = parse_date(e[4])
        doc['metadata']['sample']['host'] = e[5]
        doc['metadata']['misc']['host_gb'] = e[6]
        doc['metadata']['sample']['country'] = e[7]
        doc['metadata']['misc']['moltype'] = e[8]

        # print(json.dumps(doc, sort_keys=True, indent=4))
        try:
            db.coronaviridae.insert_one(doc)
        except pymongo.errors.DuplicateKeyError:
            continue
        bar.update(counter)


collection = db.get_collection('coronaviridae')
