'''
zoo command line.
'''

import click
from .cli.cell import init, add, commit, diff, pull, status, validate
from .cli.cell import drop, destroy
from .cli.minhash import minhash, sbt_index
from .cli.io import load, dump
from .cli.digest import digest
from .cli.scale import scale


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
cli.add_command(validate)  # TODO

# minhash/ SBT related
cli.add_command(minhash)
cli.add_command(sbt_index)

# io
cli.add_command(load)  # TODO
cli.add_command(dump)  # TODO

# alignment, tree
cli.add_command(digest)  # TODO

# port
cli.add_command(scale)  # TODO
