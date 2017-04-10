'''2 link examples.'''


import json
from zoo import get_data
from zoo.link import Link


'''
1.

We are looking at a document that contains data on the plant pest "plum
pox virus", that amongst other things infects peaches. In the "relative"
section of the document there is a link to the chloroplast sequence of the
peach reference genome. This sequence would have been too large to embed
in the document. However, the link allows random access to the indicated
fasta file, which is still fast.

Our document looks like this:

{
    "_id": "someuuid",
    "seq": "AAAATATAAA..."
    },
    "relative": {
        "link": {
            "target": "plum/plum.fa",
            "key": "NC_014697.1",
            "txt": "chloroplast prunus persica"
        }
    }
}
'''


with open(get_data('plum/plum_pox.json'), 'r+') as file:
    rec = json.load(file)

l = Link(**rec['relative']['link'])
l.is_valid()  # check if link has correct schema

# get first 10 bases of chloroplast sequence
seq = l.access_fasta(func=lambda x: x.split(' ')[0])[:10]
# "func" is what gets passed to pyfaidx.Fasta as "key_function"
assert seq == 'TGGGCGAACG'


'''
2.

Now something a bit more difficult: We'll import a seqmented virus and retieve
all segments. One way to achieve this is through a linked list, i.e. each
segment points to the next one.
'''

'''
How to represent a segmented virus in the document database? This notebook
is part tutorial, part reference.
'''


from Bio import Entrez, SeqIO
import json
from pymongo import MongoClient
from uuid import uuid4
'''
In other words, only after generating 1 billion UUIDs every second for the next
100 years, the probability of creating just one duplicate would be about 50%
(stackoverflow, 1155008).
'''


# http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc110
Entrez.email = ''  # better to provide none than a random one


'''
[article](http://www.pnas.org/content/111/18/6744.full):

Qin, X.-C. et al. A tick-borne segmented RNA virus contains genome segments
derived from unsegmented viral ancestors. PNAS 111, 6744–6749 (2014).
'''
segments = ["KJ001579", "KJ001580", "KJ001581", "KJ001582"]


l = []
for s in segments:
    handle = Entrez.efetch(
        db="nucleotide",
        id=s,  # list of GenBank identifiers
        rettype="fasta",  # rettype="gb",
        retmode="text"
        )
    l.append(SeqIO.read(handle, "fasta"))  # "genbank"


# m = []
# for s in segments:
#     handle = Entrez.efetch(
#         db="nucleotide",
#         id=s,  # list of GenBank identifiers
#         rettype="gb",  # rettype="gb",
#         retmode="text"
#         )
#     m.append(SeqIO.read(handle, "genbank"))


handle.close()


# stackoverflow, 23333123
s = m[1]
for feature in s.features:
    if feature.type == "CDS":
        print(feature.location.extract(s).seq[:20])


# http://biopython.org/wiki/SeqIO
with open("jmtv.fa", "w") as output_handle:
    SeqIO.write(l, output_handle, "fasta")


schema = '''
{
    "_id": "",
    "sequence": "",
    "link": ""
}
'''


d = json.loads(schema)


client = MongoClient("localhost:27017")
db = client["jmtv"]


# populate db
link = ""
for i in l:
    d = json.loads(schema)

    d["_id"] = str(uuid4())[:8]
    d["sequence"] = str(i.seq)[:20]
    d["link"] = link

    db.jmtv.insert_one(d)
    link = d["_id"]


collection = db.get_collection('jmtv')
collection.find_one()


a = db.jmtv.aggregate([
    {"$match": {"_id": "87847a5b"}},
    {"$graphLookup":
        {
            "from": "jmtv",
            "startWith": "$link",
            "connectFromField": "link",
            "connectToField": "_id",
            # "maxDepth": 2,  # defaults to all, first step is 0
            "as": "segments"
        }}
    ])


next(a)['segments']
'''
[{'_id': '67d0c7a0', 'link': '', 'sequence': 'GTTAAAAGAGGCCGCCCTTT'},
 {'_id': '88b696ca', 'link': '67d0c7a0', 'sequence': 'GTTTAAAAAGCGGACCGTGC'},
 {'_id': '8ac49d0b', 'link': '88b696ca', 'sequence': 'GTTAAAAAGCGCCAGCTGAG'}]
'''

collection.update_one(
    {'_id': '67d0c7a0'},
    {'$set': {'link': '87847a5b'}}
    )

'''
Create a circular linking.

In [118]: [print(i) for i in a]
{'link': '87847a5b', 'sequence': 'GTTAAAAGAGGCCGCCCTTT', '_id': '67d0c7a0'}
{'link': '67d0c7a0', 'sequence': 'GTTTAAAAAGCGGACCGTGC', '_id': '88b696ca'}
{'link': '88b696ca', 'sequence': 'GTTAAAAAGCGCCAGCTGAG', '_id': '8ac49d0b'}
{'link': '8ac49d0b', 'sequence': 'GTTTAAAAACGGCCGCCCTG', '_id': '87847a5b'}
'''
a = db.jmtv.aggregate([
    {"$match": {"_id": "87847a5b"}},
    {"$graphLookup":
        {
            "from": "jmtv",
            "startWith": "$link",
            "connectFromField": "link",
            "connectToField": "_id",
            # "maxDepth": 2,  # defaults to all, first step is 0
            "as": "segments"
        }}
    ])
next(a)['segments']
'''
[{'_id': '87847a5b', 'link': '8ac49d0b', 'sequence': 'GTTTAAAAACGGCCGCCCTG'},
 {'_id': '67d0c7a0', 'link': '87847a5b', 'sequence': 'GTTAAAAGAGGCCGCCCTTT'},
 {'_id': '88b696ca', 'link': '67d0c7a0', 'sequence': 'GTTTAAAAAGCGGACCGTGC'},
 {'_id': '8ac49d0b', 'link': '88b696ca', 'sequence': 'GTTAAAAAGCGCCAGCTGAG'}]

So it doesn't loop through graph forever, good.
'''



