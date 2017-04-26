'''
zoo command line.
'''

import click
from .cli.cell import init, add, commit, diff, pull, status
from .cli.cell import drop, destroy
from .cli.minhash import minhash, sbt_index
from .cli.io import load, dump
from .cli.digest import digest
# from .cli.sample import sample, query
from .cli.scale import scale
from .cli.schema import schema


@click.group()
def cli():
    pass


# general cell handling
cli.add_command(init)
cli.add_command(add)
cli.add_command(commit)
cli.add_command(diff)
cli.add_command(pull)
cli.add_command(status)
cli.add_command(drop)
cli.add_command(destroy)

# minhash/ SBT related
cli.add_command(minhash)
cli.add_command(sbt_index)

# io
cli.add_command(load)
cli.add_command(dump)

# schema
# port
cli.add_command(schema)  # TODO

# alignment, tree
cli.add_command(digest)  # TODO

# port
cli.add_command(scale)  # TODO

# sample, query?
# cli.add_command(sample)  # TODO
# cli.add_command(query)  # TODO


