'''
Turn into script:

1. Download from ftp.
2. Run.
3. Save database as shippable JSON.
'''

from copy import deepcopy
# from itertools import islice
import json
import pandas as pd
import progressbar
from pyfaidx import Fasta
from pymongo import MongoClient
import re
from uuid import uuid4
from zoo import get_schema
from zoo.interval import chunks
from zoo.parse import parse_date, parse_nomenclature_iav, parse_location
from zoo.utils import deep_get, deep_set


with open(get_schema('influenza_a_virus.json')) as infile:
    schema = json.load(infile)


client = MongoClient('localhost:27017')
db = client['mock']


'''
1. Load metadata.
'''


fp = 'genomeset.dat'
# see README in fp for header description
header = (
    "genbank",
    "host",
    "segment_number",
    "subtype",
    "country",
    "date",
    "seqlen",
    "isolate",
    "age",
    "gender",
    "id"
    )
df = pd.DataFrame.from_csv(fp, sep='\t', header=None, index_col=None)
df.columns = header

# fp_geo = '/Users/pi/data/geolocation/allCountries.redux.txt'
# geo = pd.DataFrame.from_csv(fp_geo, sep='\t', header=None, index_col=None)
# geo.columns = 'name lat long'.split(' ')

'''
We'll create the following document using a DotMap object. Note however that
if we wanted to enforce the schema, this is not ideal, because DotMap allows
assignment to non-existing keys. If schema enforcment is desired, look at
zoo.utils.deep_set() and deep_get(). They have been expicitly built for this
use case.
'''


bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
counter = 0

# for i in islice(df.iterrows(), 50):
for i in df.iterrows():
    counter += 1
    j = i[1]
    if 'influenza b' not in j.isolate.lower():

        d = deepcopy(schema)  # important, otherwise schema is modified

        date = parse_date(j.date)
        # try:  # date can be NaN
        #     date = parse_date(j.date)
        # except AttributeError:
        #     date = {'d': None, 'm': None, 'y': None}
        try:
            host = j.host.lower()
        except AttributeError:
            host = None

        entries = {
            '_id': str(uuid4()),
            'metadata.location': j.country,
            'metadata.date': date,
            'metadata.host': host
            }

        for k, v in entries.items():
            try:  # NaN in host, date
                deep_set(d, k, v, replace=True)
            except AttributeError:
                pass

        deep_get(d, 'metadata.alt_id').append({'genbank': j.genbank})
        deep_get(d, 'metadata.grp_id').append({'segments': j.id})

        # try:  # host: NaN
        #     deep_set(d, 'metadata.host', j.host.lower(),)
        # except AttributeError:
        #     pass

        deep_set(d, 'relative.taxonomy.subtype', j.subtype, force=True)
        deep_set(d, 'derivative.segment_number', j.segment_number, force=True)
        deep_set(d, 'derivative.length', j.seqlen, force=True)
        deep_set(d, 'metadata.age', j.age, force=True)
        deep_set(d, 'metadata.gender', j.gender, force=True)
        deep_set(
            d, 'relative.taxonomy.nomenclature',
            re.search('\((.*)\)', j.isolate).group(1))
        # format: 'Influenza A virus (A/Hong Kong/1/1968(H3N2))'
        # returns: 'A/Hong Kong/1/1968(H3N2)'
        # stackoverflow, 15864800

        # nomenclature entries such as 'A/X-31(H3N2)'
        try:
            deep_set(
                d,
                'metadata.host_detail',
                parse_nomenclature_iav(
                    deep_get(
                        d, 'relative.taxonomy.nomenclature'
                        ))['host'].lower(),
                force=True)
        except AttributeError:
            deep_set(d, 'metadata.host_detail', None, force=True)

        # insert into DB "mock" collection "iav"
        db.iav.insert_one(d)
        bar.update(counter)

# db.iav.find_one()
# db.iav.count()
# print(json.dumps(d, indent=4))


# Size increase compared to raw data: from 30 (uncompressed csv) to 50 mb
# (compressed BSON). Not so frugal.


'''
2. Explore some.
'''


db.iav.distinct('relative.taxonomy.subtype')


'''
3. Load annotation.
'''


fp = '/Users/pi/data/virus/influenza_a/raw/ncbi/influenza.dat'

with open(get_schema('annotation.json')) as infile:
    anno = json.load(infile)

# schmema nor deprecated
# {
#     "end": null,
#     "fuzzy": null,
#     "id": null,
#     "syn": "",
#     "name": "",
#     "source": "",
#     "start": null
# }

# anno = '''
# {
#     "end": null,
#     "fuzzy": null,
#     "id": null,
#     "name": "",
#     "source": "",
#     "start": null
# }
# '''


# First, fill the annotation schema.
d_anno = {}
with open(fp, 'r+') as infile:
    counter = 0  # in total 8 weird entries for location: "<1"
    for line in infile:

        entry = line.strip().split('\t')
        acc_nt = entry[0]

        l = []
        for acc_pr, coord in chunks(entry[1:], 2):
            d = deepcopy(anno)
            c = coord.split(':')[1].strip(')').split(', ')[0]
            d['source'] = 'genbank'
            # d['location'] = c  # redundant, only here for reference

            try:
                d['start'], d['end'], d['fuzzy'] = parse_location(c)
            except ValueError:
                counter += 1
                continue

            d['id'] = acc_pr
            # this is referred to as "fuzzy location"
            l.append(d)

        d_anno[acc_nt] = l


'''
...
'CY102078': [{'end': '707',
   'fuzzy': None,
   'id': 'AET75802',
   'location': '15-707',
   'name': '',
   'source': 'genbank',
   'start': '15'},
  {'end': '44',
   'fuzzy': None,
   'id': 'AET75803',
   'location': '15-44',
   'name': '',
   'source': 'genbank',
   'start': '15'}],
...
'''


# The following yould be dog slow w/o an index.
# Index alt_id genbank field, bc/ we'll use it lots.
db.iav.create_index('metadata.alt_id.genbank')
# next(db.iav.find({'metadata.alt_id.genbank': 'AJ404630'}))


# Add the annotation schema to existing records.
bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
counter = 0
for k, v in d_anno.items():
    counter += 1
    for i in v:
        db.iav.update_one(
            {'metadata.alt_id.genbank': k},
            {'$push': {'derivative.annotation': i}}
            )
    bar.update(counter)
# print(json.dumps(db.iav.find_one(), indent=4))


# We'll now insert the corresponding name to a given annotation ID.
# But first, we'll index the ID field.
db.iav.create_index('derivative.annotation.id')


# Import annotation names.
fp = '/Users/pi/data/virus/influenza_a/raw/ncbi/influenza_aa.dat'
df = pd.DataFrame.from_csv(
    fp, sep='\t', header=None, index_col=None)[[0, 2]]
df.columns = ['id', 'name']


bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
counter = 0
for i in df.iterrows():
    counter += 1
    k, v = i[1]
    db.iav.update_one(
        {'derivative.annotation.id': k},
        {'$set': {'derivative.annotation.$.name': v}}
        # problem: set does not work on array, stackoverflow, 23821392
        # illusion to work for couple of entries, but those are Influenza B
        # include .$.
    )
    bar.update(counter)
# print(json.dumps(db.iav.find_one(), indent=4))


'''
4. Load sequences.
'''

fp = '/Users/pi/data/virus/influenza_a/raw/ncbi/influenza.fna'
fa = Fasta(fp, key_function=lambda x: x.split('|')[3])
# key w/o key_function: 'gi|59296|gb|X66320|Influenza'
# gen = (i for i in fa.keys())
# k = next(gen)
# 'X66320'


bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
counter = 0
for i in fa:
    counter += 1
    db.iav.update_one(
        {'metadata.alt_id.genbank': i.name},  # already indexed
        {'$set': {'sequence': str(i)}}
    )
    bar.update(counter)


'''
GISAID data, include. Do string match or partial string match (if GISAID only
saved coding sequences) to integrate the two resources.

Complex query: Get all non-coding sequences.

Develop EDA functions.
'''

























