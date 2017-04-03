'''
get a number of docs, to SBT in tempfile or disk, query an interesting
seq, save result in doc?

sbt_index([cursor], fp='', temporary=True)
sbt_combine([list, of, cursors], fp='', temporary=True)
sbt_search(sbt, query)
'''


# https://github.com/luizirber/2016-sbt-minhash/blob/master/notebooks/SBT%20with%20MinHash%20leaves.ipynb
# http://blog.luizirber.org/2016/12/28/soursigs-arch-1/
from pyfaidx import Fasta
from sourmash_lib import Estimators
from sourmash_lib.sbt import SBT, GraphFactory
from sourmash_lib.sbtmh import SigLeaf, search_minhashes
from sourmash_lib.signature import SourmashSignature
from zoo import get_data


KSIZE = 16
N = 200

with open(get_data('data/mock_flu.fa'), 'r+') as file:
    fa = Fasta(file)


# init SBT
factory = GraphFactory(ksize=KSIZE, starting_size=1e5, n_tables=4)
# 4 .. nt?
tree = SBT(factory, d=2)  # d .. see "n-ary " in notebook
'''
How does the internal nodes / total ratio affect query times? Test and
put as section in master thesis.
'''


# load signatures to tree
for i in fa:
    key = i.name
    seq = str(fa[key])  # db.ref.find_one()['sequence']  # 'ACTG...'
    e = Estimators(ksize=KSIZE, n=N)
    e.add_sequence(seq)  # e.get_hashes()
    s = SourmashSignature(email='', estimator=e, name=key)
    # s.estimator.get_hashes()
    # s.name()
    leaf = SigLeaf(metadata=key, data=s)
    # SigLeaf(metadata, data, name=None)
    tree.add_node(node=leaf)
# tree.print()  # ignore pylint


# search the last fasta entry against the SBT (">0.95")
filtered = tree.find(search_minhashes, s, 0.1)
matches = [(str(i.metadata), i.data.similarity(s)) for i in filtered]
# [('0.95', 1.0)]  # fasta header, similarity


tree.save('mock_flu')


'''shell
head -n2 mock_flu.fa | sourmash compute -k 16 -n 200 -o virion.json -
sourmash sbt_search -k 16 mock_flu.sbt.json virion.json
# header: similarity, fasta header ("key" above)
# 1.00 0.00
# 0.44 0.05
# 0.34 0.10
# 0.14 0.15
# 0.10 0.20
'''









