'''
- load core schema
- load entries from MSA and validate inserted files
- encode MSA and decode into original file: md5 hash to prove identity
- commit, dat share
- EDA w/ metadata: map, time series
- use annotations: number of mutations per site in e.g. RdRp
- encode tree (modification)
- reference free variant calling and storing result in relative
- dat pull, diff, add, status
- reconstruct and draw tree
- create registry log w/ zoo status, search the minhash signature attached to
this against a SBT made from all reference sequences
- load all other Ebola data from VIPR into another collection
- ML: use both collections for host prediction
- random access fasta experiments here as well (linked host genome)
'''


# from itertools import islice
import json
from jsonschema import validate
import pandas as pd
from pyfaidx import Fasta
from pymongo import MongoClient
from uuid import uuid4
from zoo import get_schema, get_data
from zoo.align import encode_gaps, decode_gaps, hash_seq
from zoo.utils import deep_set, deep_get


client = MongoClient('localhost:27017')
c = client['ebola']['makona']


fp = get_data('ebola/Makona_1610_genomes_2016-06-23.fasta')
fp_geo = get_data('ebola/geo_lat_long.tsv')
fp_makona = get_data('ebola/Makona_1610_metadata_2016-06-23.csv')


'''
load core schema
'''


with open(get_schema('core.json'), 'r') as file:
    core = json.load(file)
with open(get_schema('metadata.json'), 'r') as file:
    metadata = json.load(file)
with open(get_schema('relative.json'), 'r') as file:
    relative = json.load(file)


deep_set(core, 'properties.meta', metadata, replace=True)
deep_set(core, 'properties.rel', relative, replace=True)


'''
load entries from MSA and validate documents before insertion
'''


fa = Fasta(fp, key_function=lambda x: x.split('|')[1])
h = hash_seq(str(i) for i in fa)


# get location names right
df_geo = pd.read_csv(fp_geo, sep='\t', header=0)
df_makona = pd.read_csv(fp_makona, sep=',', header=0)

s = set()
for i in df_makona['location']:
    try:
        s.update([i.lower()])
    except AttributeError:  # 'float' object has no attribute 'lower' for nan
        pass


# s - set(df_geo['location'])  # manually determine synonyms
syn = {
    '?': '',
    'grandbassa': 'grand_bassa',
    'grandcapemount': 'grand_cape_mount',
    'grandkru': 'grand_kru',
    'portloko': 'port_loko',
    'rivergee': 'river_gee',
    'westernarea': 'western_area',
    'westernrural': 'western_rural',
    'westernurban': 'western_urban',
    'SLE': 'sierra_leone',
    'LBR': 'liberia',
    'GIN': 'guinea'
    }

# for i in islice(fa, 3):
for i in fa:
    _, strain, someid, country, location, date = i.long_name.split('|')
    seq = str(i)

    # consistent country names to assign long, lat
    country = syn[country]
    location = location.lower()
    if location == '?':
        location = country
    elif location in syn:
        location = syn[location]

    try:
        lat, long = df_geo[
            df_geo['location'] == location
            ][['latitude', 'longitude']].values[0]
        lat = round(float(lat), 3)
        long = round(float(long), 3)
    except IndexError:
        lat, long = None, None

    rec = {}  # record
    rec['_id'] = str(uuid4())

    # metadata
    deep_set(rec, 'meta.geo.cou', country, force=True)
    deep_set(rec, 'meta.geo.loc', location, force=True)
    deep_set(rec, 'meta.geo.lat', lat, force=True)
    deep_set(rec, 'meta.geo.long', long, force=True)
    deep_set(rec, 'meta.ids', [], force=True)
    deep_get(rec, 'meta.ids').append({'sample': someid})

    # y, m, d = map(int, date.split('-'))
    # ymd = {k: v for k, v in zip('ymd', map(int, date.split('-')))}
    deep_set(rec, 'meta.date', date, force=True)

    # relative
    deep_set(rec, 'rel.tax.strain', strain, force=True)
    deep_set(rec, 'rel.tax.species', 'ebola', force=True)

    # more metadata
    # - Makona_1610_metadata_2016-06-23.csv

    '''
    encode MSA
    '''

    # gaps and MSA hash as ID
    seq_true = seq.replace('-', '')
    gaps = encode_gaps(seq)
    # ignores row order

    rec['seq'] = seq_true
    deep_set(rec, 'rel.msa', [], force=True)
    deep_get(rec, 'rel.msa').append(
        {'hash': h, 'gaps': gaps, 'txt': 'Makona_1610_genomes_2016-06-23'}
        )

    # check record against schema
    assert validate(rec, core) is None

    c.insert_one(rec)
c.count()  # 1610


# All location entries should be present in geolocation df:
set(c.distinct('meta.geo.loc')) - set(df_geo['location'])  # {''}


'''
decode MSA, md5 hash should be equal before/ after en-/ decoding
'''


q = c.find(
    {'rel.msa.hash': 'a130dfeee97e58bf04079d9efe358a8b'},
    {'seq': 1, '_id': 0, 'rel.msa.gaps': 1})


gen = (decode_gaps(i['seq'], deep_get(i, 'rel.msa')[0]['gaps']) for i in q)
print(hash_seq(gen))  # a130dfeee97e58bf04079d9efe358a8b





















