'''2 link examples.'''


from Bio import Entrez, SeqIO
from copy import deepcopy
import json
from pymongo import MongoClient
from zoo import get_data
from zoo.link import Link
from zoo.utils import deep_get


'''
1.

We are looking at a document that contains data on the plant pest 'plum
pox virus', that amongst other things infects peaches. In the 'relative'
section of the document there is a link to the chloroplast sequence of the
peach reference genome. This sequence would have been too large to embed
in the document. However, the link allows random access to the indicated
fasta file, which is still fast.

Our document looks like this:

{
    '_id': 'someuuid',
    'seq': 'AAAATATAAA...'
    },
    'relative': {
        'link': {
            'target': 'plum/plum.fa',
            'key': 'NC_014697.1',
            'txt': 'chloroplast prunus persica'
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
# 'func' is what gets passed to pyfaidx.Fasta as 'key_function'
assert seq == 'TGGGCGAACG'


'''
2.

Now something a bit more difficult: We'll import a seqmented virus and retieve
all segments. One way to achieve this is through a linked list, i.e. each
segment points to the next one.

How to represent a segmented virus in the document database? This notebook
is part tutorial, part reference.

[article](http://www.pnas.org/content/111/18/6744.full):

Qin, X.-C. et al. A tick-borne segmented RNA virus contains genome segments
derived from unsegmented viral ancestors. PNAS 111, 6744–6749 (2014).
'''
segments = ['KJ001579', 'KJ001580', 'KJ001581', 'KJ001582']


l = []
Entrez.email = ''
for s in segments:
    handle = Entrez.efetch(
        db='nucleotide',
        id=s,  # list of GenBank identifiers
        rettype='fasta',  # rettype='gb',
        retmode='text'
        )
    l.append(SeqIO.read(handle, 'fasta'))
handle.close()


with open('jmtv.fa', 'w+') as output_handle:
    SeqIO.write(l, output_handle, 'fasta')
# This file can be found in zoo/data/.


# document schema
schema = json.loads('''
    {
        "_id": "",
        "seq": "",
        "rel": {
            "link": null
        }
    }
    ''')


c = MongoClient('localhost:27017')['test']['jmtv']
key = ''
for i in l:
    rec = deepcopy(schema)
    rec['_id'] = i.name
    rec['seq'] = str(i.seq)[:10]  # for brevity, we'll use only 10 bases
    rec['rel']['link'] = []
    rec['rel']['link'].append(
        Link(target='test.jmtv._id', key=key).__dict__)
    '''
    We don't want to pass the object itself, but the dict representation.
    This is more convenient, when sending records around, bc/ dicts are
    trivial to translate to JSON.
    '''
    if Link(**rec['rel']['link'][0]).is_valid():  # ugly
        c.insert_one(rec)
    key = i.name  # pass the previous segments ID as link key to the next seq


'''
One nice way to store a segmented genome is by a circular linked list where the
last element connects to to first. That way, we do not store too many links
(e.g. compared to saving all n-1 segment links in each of the n documents).
MongoDB with its graph database functionality is intelligent enough to not
cycle indefiniely though such a circular list, but stops once all nodes have
been visited, i.e. once all segments have been collected.
'''
# create circular linked list
c.update_one(
    {'_id': 'KJ001579.1'},
    {'$set': {'rel.link': [{'key': 'KJ001582.1', 'target': 'test.jmtv._id'}]}}
    )


c.find_one()
'''
{'_id': 'KJ001579.1',
 'rel': {'link': [{'key': 'KJ001582.1', 'target': 'test.jmtv._id'}]},
 'seq': 'GTTAAAAGAG'}
'''
c.count()
# 4


'''
We'll now do two things:

1.1 use MongoDB's "$graphLookup" operator to travers the "graph" (linked list)
1.2 access a record linked from another record (zoo functionality)
'''

'''
1.1 Perform a graph traversal to get aggregate all linked elements.
'''

a = c.aggregate([
    {'$match': {'_id': 'KJ001579.1'}},
    {'$graphLookup':
        {
            'from': 'jmtv',  # make sure this is really the collection
            'startWith': '$rel.link.key',
            'connectFromField': 'rel.link.key',
            'connectToField': '_id',
            'maxDepth': 100,  # defaults to all, first step is 0
            'as': 'segments'
        }}
    ])

print([i['_id'] for i in next(a)['segments']])
'''
# ['KJ001580.1', 'KJ001579.1', 'KJ001581.1']

Note how even though we allow a maximum recursion depth of 100, the graph
traversal stops after all segments are recovered.
'''

'''
1.2 Access individually linked document.
'''

rec = next(c.find({'_id': 'KJ001580.1'}))
'''
link
# key: KJ001579.1
# target: test.jmtv._id
'''
link = Link(**deep_get(rec, 'rel.link')[0])
next(link.access_document())
'''
{'_id': 'KJ001579.1',
 'rel': {'link': [{'key': 'KJ001582.1', 'target': 'test.jmtv._id'}]},
 'seq': 'GTTAAAAGAG'}
'''
