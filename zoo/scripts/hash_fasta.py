'''Hash fasta sequences.

A fasta file is transformed into a fasta file with the header containing (only)
the hash (md5, hardcoded) of a given entry's sequence. Another file (yaml
format) is created containing the mapping from the original header to the
hash.

The intended use is to fill a (document-based NoSQL) database with the metadata
in the fasta header but to leave the sequences in flat files. We can then
query the database to return hash ids of interest, and use a fasta index (
samtools, pyfaidx) to access the sequences from the (hashed, indexed) fasta
file. This advantage are manyfold:

* we have a universal id, which does not depend on the fact that a sequence
  has een assigned an accession id, which might or might not be consistent
  accross databases (e.g. for influenza, there is no one-to-one mapping for
  sequence data from GenBank and GISAID)
* duplicate sequences are easily identified (and removed)
* the database is not cluttered with sequence data (for viruses nucleotide
  sequences have a trivial size, but bacteria already pose a problem to
  scalability)
* duplicate sequences get removed from flat fasta
* the fasta files can be structured independently from the database, e.g.
  in distributed blocks accross a cluster, which is scalable

Example:
  python hash_fasta.py in.fa out
  # Duplicates?
  grep ">" test.fa | wc -l
  # 163080
  grep ">" test.fa | sort | uniq | wc -l
  # 163080
  # No.

Usage:
  hash_fasta.py <infile> <outprefix>
  hash_fasta.py (-h | --help)
  hash_fasta.py --version

Arguments:
    infile                      Fasta file with arbitrary header.

Options:
  -h --help                     Show this screen.
  --version                     Show version.
'''

from docopt import docopt
import hashlib
import progressbar
import skbio


if __name__ == '__main__':
    arguments = docopt(
        __doc__, version='0.1dev')


infile = arguments['<infile>']
outprefix = arguments['<outprefix>']
outfasta = outprefix + '.fa'
outyaml = outprefix + '.yaml'

x = skbio.io.read(infile, format='fasta')


duplicates = set()
bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
counter = 0
with open(outfasta, 'w+') as outf, open(outyaml, 'w+') as outm:

    print('Hashing things:')
    for i in x:
        s = str(i)
        se = s.encode('ASCII')
        m = hashlib.md5()  # sha256(), ...
        m.update(se)
        h = m.hexdigest()

        outm.write(i.metadata['id'] + ': ' + h + '\n')
        # header info: 6409dac0c80c79ddf61858be7f2c699b

        if h not in duplicates:
            outf.write('>' + h + '\n' + s + '\n')
            # >f907632441d101ec5a84150594a88a1b
            # ACTGA...
            duplicates.update([h])

        counter += 1
        if counter % 100 == 0:
            bar.update(counter)

print('\nDone.')
