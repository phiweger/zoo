import click
from sourmash_lib import Estimators
from sourmash_lib.sbt import SBT, GraphFactory
from sourmash_lib.sbtmh import SigLeaf
from sourmash_lib.signature import SourmashSignature


@click.command()
def sbt_index():
    '''Create a sequence Bloom tree from a cell/ database cursor.
    1. select seqs for tree
    2. assign common id (field derivative.minhash.sbt.ids)
    3. minhash seqs, name == UUID, md5? (think about SBT reuse)
    4. query a different collection/ metagenome against this

    --index {raw, minhash}
    input: all of cell or cursor
    '''
    pass


@click.command()
def minhash():
    '''Minhash a cell/ database cursor.
    just plain old sigs for collection
    '''
    pass
