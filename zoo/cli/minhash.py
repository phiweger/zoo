import click
from sourmash_lib import Estimators
from sourmash_lib.sbt import SBT, GraphFactory
from sourmash_lib.sbtmh import SigLeaf
from sourmash_lib.signature import SourmashSignature


@click.command()
def sbt_index():  # drop database entirely
    print('Trying.')
