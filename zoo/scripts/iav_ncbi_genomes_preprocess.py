
'''
To understand why NCBI Genomes deposits three files, think of how you would
model this data in a relational database. They then take this scheme and dump
it awkwardly, so we have to put it back together.

acc_nt .. GenBank accession ID for nucleotide sequence
acc_pr .. GenBank accession ID for protein sequence
coord .. coordinates

- map_ntpr .. maps accession nt to (multiple) accessions pr,
including annotation coords
- map_pr .. pr accession, pr name
- map_nt .. nt accession, nt name

1. Create a (raw) document structure. Data need not be cleaned yet.
2. Filter:
    - date
    - location
    - influenza name (exclude IBV)
    - extract host and specific host (avian, duck), all lowercase
    - use fauna (nextflu) country to get geodata
3. Load to mongodb.
4. Explore.
'''


# from collections import Counter
# from itertools import islice
# from dotmap import DotMap
# from pprint import pprint
from pyfaidx import Fasta
from pymongo import MongoClient
import progressbar
import sys
sys.path.append("/Users/pi/repositories/")
from toolshed.utils_date import parse_date


def chunks(l, n):
    '''
    Yield successive n-sized chunks from l (stackoverflow, 312443).

    a = [1, 2, 3, 4]
    list(chunks(a, 2))
    # [[1, 2], [3, 4]]

    Returns empty list if list empty.
    '''
    for i in range(0, len(l), n):
        yield l[i:i + n]


def locate(locations, seq):
    '''
    Return the (concatenated) sequence string corresponding to the
    list of locations (a location being a tuple of start and end
    coordinates).

    sample.annotation.M1.location
    # ['6-764']

    764 - 6 + 1
    # 759

    # len(locate(sample.annotation.M1.location, sample.sequence))
    759

    # checked that it works for 2 locations as well, should be do so for > 2
    '''
    import re

    product = ''
    for location in locations:
        start, end = re.findall('\\d+', location)
        fragment = seq[int(start) - 1:int(end)]
        # -1 because GenBank entries start indexing at 1, not 0 as in python
        # when slicing e.g. [1:5], python excludes position 5, so here no -1
        product += fragment
        # order matters, eyeballed some entries and they seem ordered
    return product


fp_nt = '/Users/pi/data/influenza/raw/ncbi_genomes/genomeset.dat'
fp_pr = '/Users/pi/data/influenza/raw/ncbi_genomes/influenza_aa.dat'
fp_ntpr = '/Users/pi/data/influenza/raw/ncbi_genomes/influenza.dat'
fp_fa = '/Users/pi/data/influenza/raw/ncbi_genomes/influenza.fna'


# see README in fp for header description
header = (
    "accession_genbank",
    "host",
    "segment_number",
    "subtype",
    "country",
    "ymd",
    "seqlength",
    "isolate",
    "age",
    "gender",
    "id"
    )


'''
0. TODO

- parse_iav function
- apply annotation to sequence (use Biopython?)
'''


'''
1. Create a (raw) document structure. Data need not be cleaned yet.
'''


map_nt = {}
with open(fp_nt, 'r+') as nt:
    # for line in islice(genomeset, 10):
    for line in nt:
        if line != '\n':  # weird file format w/ \n between isolate entries
            d = {}
            entry = line.strip().split('\t')
            z = zip(
                header[1:],
                entry[1:]
                )
            for key, value in z:
                d[key] = value
            map_nt[entry[0]] = d


map_pr = {}
with open(fp_pr, 'r+') as pr:
    for line in pr:
        entry = line.strip().split('\t')
        acc_pr, name = [entry[i] for i in (0, 2)]
        map_pr[acc_pr] = name

# KY171142
map_ntpr = {}
with open(fp_ntpr, 'r+') as ntpr:
    for line in ntpr:

        entry = line.strip().split('\t')
        acc_nt = entry[0]

        l = []
        for acc_pr, coord in chunks(entry[1:], 2):
            name = map_pr[acc_pr]
            d = {}

            c = coord.split(':')[1].strip(')').split(', ')
            d['name'] = name
            d['accession'] = acc_pr
            d['location'] = c
            d['fuzzy'] = int(any([i in coord for i in ['>', '<']]))
            # this is referred to as "fuzzy location"
            l.append(d)

        map_ntpr[acc_nt] = l


fa = Fasta(fp_fa, key_function=lambda x: x.split('|')[3])
# imports original header 'gi|90572149|gb|CY010182|Influenza' as 'CY010182'


documents = {}
for key in map_nt.keys():
    documents[key] = {}
    documents[key]['annotation'] = map_ntpr[key]
    documents[key]['metadata'] = map_nt[key]
    documents[key]['sequence'] = str(fa[key])
    documents[key]['_id'] = key  
    # index of a document in mongodb
    # _id means GenBank (nucleotide) accession number
    # the accession number in annotation refers to the protein


'''
2. Filter:
    - date
    - location
    - influenza name (exclude IBV)
    - extract host and specific host (avian, duck), all lowercase
    - use fauna (nextflu) country to get geodata

how nextflu [does it](https://github.com/blab/nextflu/tree/master/augur)

> Filter: Keeps viruses with fully specified dates, cell passage and only one
sequence per strain name. Subsamples to 50 (by default) sequences per month for
the last 3 (by default) years before present. Appends geographic metadata.
Subsampling prefers longer sequences over shorter sequences and prefer more
geographic diversity over less geographic diversity.
'''


# leave for later, move on to 3. for now


'''
3. Load to mongodb.
'''


# from pymongo import MongoClient
# import progressbar
client = MongoClient("localhost:27017")
db = client["zoo"]


'''
# in mongodb shell
use zoo
show dbs
show collections
db.influenza_a_virus.find()  # display 20 results
db.influenza_a_virus.findOne()  # pretty print
db.runCommand( { dropDatabase: 1 } )  # delete database

# look at "schema"/ data model of a collection
# stackoverflow, 14713179
db.collectionName.find().pretty()

# backup (binary)
# stackoverflow, 4880874
mongodump --db zoo --collection influenza_a_virus --out - | \
gzip > collectiondump.gz
# or w/ date
mongodump --db zoo --collection influenza_a_virus --out - | \
gzip > dump_`date "+%Y-%m-%d"`.gz
# 130 mb

# backup (valid JSON)
mongoexport -d zoo -c influenza_a_virus -o tmp/influenza_a_virus.json \
--jsonArray --pretty
# 530 mb
'''


bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
count = 0
for key, document in documents.items():
    count += 1
    if 'Influenza A virus' in document['metadata']['isolate']:
        db.influenza_a_virus.insert_one(document)
    else:
        db.influenza_b_virus.insert_one(document)
    bar.update(count)
# / 293648 Elapsed Time: 0:03:33


'''
Update.
'''


query = col.find()  # select all entries
l = []
for i in query:
    l.append(i['_id'])


bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
counter = 0
for j in l:
    q = col.find({'_id': j})
    date = next(q)['metadata']['ymd']
    col.update_one({'_id': j}, {'$set': {'metadata.date': parse_date(date)}})
    counter += 1
    bar.update(counter)


'''
Delete.
'''


# col.update_many({}, {'$unset': {'date': 1}}, False, True)  # accidental
col.update_many({}, {'$unset': {'metadata.ymd': 1}}, False, True)


'''
4. Explore in mongodb.
'''

# ..., see usage.py in ~/project/zoo


# # now retieve data
# col = db.get_collection('influenza_a_virus')
# col.find_one()
# col.find({'subtype': 'H3N2'})

# # stackoverflow, 16002659
# a = col.find({'metadata.subtype': 'H1N1'})
# a.count()
# a.next()  # like generator, yields document as a dict

# # from dotmap import DotMap
# sample = DotMap(a.next())
# sample.annotation.PB2.location
# # ['28-2307']

# sample = DotMap(a.next())
# locate(sample.annotation.M1.location, sample.sequence)
# # ATG...

# rdrp = 'PB1 PB2 PA'.split(' ')
# query = col.find({
#     'metadata.subtype': 'H1N1',
#     'annotation.name': {"$in": rdrp},
#     'annotation.fuzzy': 1
#     })
# # returns very instructive example, nice illustration of db use

# # from pprint import pprint
# # query.count()
# # for i in islice(query, 1):
# #     pprint(i)

# sample = query.next()
# for i in sample['annotation']:
#     print(i['name'] + '\n' + locate(i['location'], sample.sequence)[:10])
# # PA
# # ATGGAAGACT
# # PA
# # AATGCATCCT




# # '$exists is called an operator', there are more, see docs
# # stackoverflow, 23701414
# # stackoverflow, 2495932
# # stackoverflow, 7431422

# # Does locate() work with fuzzy coordinates?
# a = col.find({
#     'metadata.subtype': 'H1N1',
#     'annotation.M2.fuzzy': 1
#     })
# sample = DotMap(next(a))
# locate(sample.annotation.M2.location, sample.sequence)
# # ATG...
# # It does.


# '''
# Example query: Find all RdRp segments in the database and retain only those
# of IAV instances where at least 8 Segments are present.
# '''

# countries = 'Germany Thailand'.split(' ')
# query = col.find({
#     'metadata.subtype': 'H1N1',
#     'metadata.country': {"$in": countries}
#     })

# # .limit(1).toArray()
# # stackoverflow, 10885044

# rdrp = 'PB1 PB2 PA'.split(' ')
# query = col.find({
#     'metadata.subtype': 'H1N1',
#     'annotation.PA': {"$exists": True}
#     })






# query = col.distinct('metadata.country')


# '''
# -------------------------------------------------------------------------------
# '''


# # l = []
# # for key, value in documents.items():
# #     l.append(documents[key]['metadata']['isolate'])


# # m = []
# # for i in l:
# #     try:
# #         m.append(i.split(' (')[1])
# #     except IndexError:
# #         continue

# # for i in m:
# #     try:
# #         parse_iav(i)
# #     except AssertionError:
# #         continue


# '''
# example entry:

# 'KM654693': {
#     'PB1': {
#         'accession': 'AIT92315',
#         'fuzzy': 0,
#         'location': ['25-2298']
#             },
#     'PB1-F2': {
#         'accession': 'AIT92316',
#         'fuzzy': 0,
#         'location': ['233-391']
#         }
#     }

# usage:

# map_ntpr['CY158482']
# {'M1': {'accession': 'AHB21854', 'fuzzy': 0, 'location': ['15-773']},
#  'M2': {'accession': 'AHB21855', 'fuzzy': 0, 'location': ['15-40', '729-996']}}
# '''


# # influenza.dat is the mapping of genbank nt accession and genbank protein acc.
# # annotation nomenclature
# # write test with  mock record with three annotations


# '''
# unique annotation entry contains: protein id, nucleotide id, coordinates, i.e.
# "On sequence <nt id> the <coords> correspond to <p id>."
# '''


# '''
# Info about GenBank coordinates:

# - 1-based
# - notation: see [here](http://www.insdc.org/files/feature_table.html#3.4.3)
# '''


# '''
# In [65]: Counter(d_protein.values())
# Out[65]:
# Counter({'5': 1,
#          '6': 2,
#          '7': 3,
#          '8': 1,
#          'BM2': 5299,
#          'CM2': 178,
#          'HA': 120090,
#          'HE': 281,
#          'M': 1,
#          'M1': 61422,
#          'M2': 54918,
#          'M42': 2,
#          'N40': 4,
#          'NA': 77145,
#          'NB': 8879,
#          'NP': 49367,
#          'NS1': 51795,
#          'NS2': 50725,
#          'NS3': 2,
#          'P3': 224,
#          'P42': 15,
#          'PA': 48362,
#          'PA-X': 21144,
#          'PB1': 48629,
#          'PB1-F2': 28828,
#          'PB2': 48756})
# '''


# '''
# Genbank seems even less curated

# from Bio import SeqIO
# # https://www.biostars.org/p/175427/


# fp = '/Users/pi/Downloads/sequence.gb'
# for gb_record in SeqIO.parse(open(fp, "r"), "genbank") :
#     # now do something with the record
#     print("Name %s, %i features" % (gb_record.name, len(gb_record.features)))
#     print(repr(gb_record.seq))

# gb_record.features
# a = gb_record.features[3]
# a.qualifiers['gene']


# fp = '/Users/pi/Downloads/sequence.gb'
# for gb_record in SeqIO.parse(open(fp, "r"), "genbank") :
#     # now do something with the record
#     print("Name %s, %i features" % (gb_record.name, len(gb_record.features)))
#     print(repr(gb_record.seq))
# '''






