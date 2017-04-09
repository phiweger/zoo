'''
Here we show a mock JSON record. Both, the record's sequence as well as
an associated genome (plum pox virus, reference sequence and prunus persica
chloroplast reference assembly) are linked (with the JSON schema fragment
of a "link"). We can nevertheless efficiently work with the underlying
sequences, because we have (random) access to the linked files.
'''


import json
from zoo import get_data
from zoo.link import link_access


with open(get_data('plum/plum_pox.json'), 'r+') as file:
    rec = json.load(file)
'''
{'_id': 'someuuid',
 'relative': {'link': {'fp': 'NC_014697.fa',
   'txt': 'chloroplast prunus persica',
   'type': 'filepath'}},
 'seq': {'link': {'fp': 'NC_001445.fa',
   'txt': 'plum pox reference genome',
   'type': 'filepath'}}}
'''


# >NC_001445.1 Plum pox virus, complete genome.
link_access(rec['seq'], 'NC_001445.1', func=lambda x: x.split(' ')[0])
# 'AAAATATAAA...'
link_access(rec['relative'], 'NC_014697.1', func=lambda x: x.split(' ')[0])
# 'TGGGCGAACG...'
