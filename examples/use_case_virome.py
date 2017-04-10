'''
Shi, M. et al. Redefining the invertebrate RNA virosphere. Nature 540, 539–543
(2016).

NCBI accession codes: KX882764 - KX884872

Aim: use RdRp anno to (try) to predict host via GBT.

From original publication:

> Data availability. All new sequence reads generated here are available at
the NCBI Sequence Read Archive (SRA) database under the BioProject accession
PRJNA318834 (Supplementary Table 1). All virus genome sequences generated in
this study have been deposited in GenBank under the accession numbers
KX882764–KX884872 (Supplementary Table 2). All viruses discovered in this study
(fasta format), sequence alignments (fasta format), and phylogenetic
trees (newick format) are available at
https://figshare.com/articles/
Redefining_the_ invertebrate_RNA_virosphere/3792972.

Aim: collect this information in one data structure and share.
'''


from Bio import Entrez, SeqIO
from progressbar import ProgressBar, UnknownLength
from pymongo import MongoClient
from urllib.error import URLError


# http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc110
Entrez.email = ''  # it is recommended to provide an address


acc = ['KX' + str(i) for i in range(882764, 884872 + 1)]
selection = acc[:100]
with open('rna_virome_shi2016.txt', 'w+') as file:
    for i in selection:
        file.write(''.join([i, '\n']))
# len 2109


# acc = acc[:5]
l = []
counter = 0
bar = ProgressBar(max_value=UnknownLength)
for i in selection:  # as a test data set
    try:
        handle = Entrez.efetch(
            db="nucleotide",
            id=i,  # list of GenBank identifiers
            rettype="gb",  # rettype="gb",
            retmode="text"
            )
    except URLError:
        print(i, 'timed out or something.')
    l.append(SeqIO.read(handle, "genbank"))  # "genbank"
    counter += 1
    bar.update(counter)


with open('test.fa', 'w+') as file:
    for i in l:
        file.write(i.format('fasta'))

record = l[0]

record.id
record.seq
record.features
f = record.features[0]
f.id
f.location
f.extract(record.seq)


# See what's annotated.
s = set()
for record in l:
    for f in record.features:
        try:
            s.add(f.qualifiers['product'][0])
        except KeyError:
            continue


# And that is one study.
syn = [
    'RNA-dependent RNA polymerase',
    'RNA-dependent RNA polymrease',
    'RdRp',
    'PA',
    'polymerase PA',
    'polymerase PB1',
    'polymerase PB2',
    'polymerase-associated protein'
]


counter = 0
s = set()
for record in l:
    for f in record.features:
        try:
            if f.qualifiers['product'][0] in syn:
                counter += 1
        except KeyError:
            continue
# counter  # 587 of 2109
           







