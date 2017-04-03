'''
Download Influenza A virus (A/WSN/1933 TS61(H1N1)) segment 6, complete
sequence. Then mutate with increasing mutation rate to generate divergent
sequences. The resulting test data set is intended to be used in the
SBT functionality tests.
'''


from Bio import Entrez, SeqIO
import numpy as np
import random
from scipy.stats import uniform


def mutate(text, error_rate, alphabet=None):
    '''
    insert changes in a text, i.e. mistakes, increasing the
    distance to the original

    alphabet=ascii_letters
    alphabet=set(text)

    don't modify alphabet, i.e. transforming the string to lower
    '''
    if alphabet is None:
        # use text'alphabet by default
        alphabet = list(set(text))
    lt = len(text)
    num_mutations = sum(uniform.rvs(loc=0, scale=1, size=lt) < error_rate)
    sample = random.sample(range(0, lt), num_mutations)
    l = list(text)
    for i in sample:
        l[i] = random.choice(alphabet)
    return(''.join(l))


Entrez.email = ''  # better to provide none than a random one
id_genbank = 'CY010790'
handle = Entrez.efetch(
    db="nucleotide",
    id=id_genbank,  # list of GenBank identifiers
    rettype="fasta",  # rettype="gb",
    retmode="text"
    )
seq = SeqIO.read(handle, "fasta")
handle.close()


mut = []
for i in np.arange(0, 1, 0.05):
    mut.append(
        (
            '{0:0.2f}'.format(i),  # stackoverflow, 22222818
            mutate(str(seq.seq), i)
            )
        )


with open('mock_flu.fa', 'w+') as fasta:
    for m in mut:  # includes the unmodified sequence bc/ np.arange includes 0
        fasta.write(
            '>{}\n{}\n'.format(*m)
            )








