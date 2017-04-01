import click
from pymongo import MongoClient
from sourmash_lib import Estimators
from sourmash_lib.sbt import SBT, GraphFactory
from sourmash_lib.sbtmh import SigLeaf
from sourmash_lib.signature import SourmashSignature


@click.option('--client', default='localhost:27017')
@click.option('--db', required=True)
@click.option('--cell', required=True, help='Cell name.')
@click.option(
    '--query', required=False, help='A JSON formatted MongoDB query.')
@click.option('--ksize', default=31)
@click.option('--nsketch', default=1000)
@click.option('--key', default='_id')
@click.argument('file', type=click.Path())
@click.command()
def sbt_index(client, db, cell, ksize, nsketch, key, out, file):
    '''Create a sequence Bloom tree from a cell/ database cursor.
    1. select seqs for tree
    2. assign common id (field derivative.minhash.sbt.ids)
    3. minhash seqs, name == UUID, md5? (think about SBT reuse)
    4. query a different collection/ metagenome against this

    --index {raw, minhash}
    input: all of cell or cursor

    TODO: add query
    '''
    # init SBT
    factory = GraphFactory(ksize=ksize, starting_size=1e5, n_tables=4)
    # 4 .. nt?
    tree = SBT(factory, d=2)  # d .. see "n-ary " in notebook

    c = MongoClient(client)[db][cell]
    for d in c.find():
        e = Estimators(ksize=ksize, n=nsketch)
        e.add_sequence(d['sequence'], force=True)
        s = SourmashSignature(estimator=e, name=d[key])
        leaf = SigLeaf(metadata=key, data=s)
        tree.add_node(node=leaf)
    tree.save(file)


@click.command()
def minhash():
    '''Minhash a cell/ database cursor.
    just plain old sigs for collection
    '''
    pass
