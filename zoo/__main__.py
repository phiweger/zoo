'''
zoo command line.
'''

import click
from .cli.cell import init, add, commit, diff, pull, status, drop, destroy, validate
from .cli.minhash import minhash, sbt_index
from .cli.io import io
from .cli.align import msa_encode, msa_decode, sam_encode, sam_decode
from .cli.tree import tree_encode, tree_decode


@click.group()
def cli():
    pass


# general cell handling
cli.add_command(init)
cli.add_command(add)
cli.add_command(commit)
cli.add_command(diff)       # TODO
cli.add_command(pull)       # TODO
cli.add_command(status)     # TODO
cli.add_command(drop)       # TODO
cli.add_command(destroy)    # TODO
cli.add_command(validate)    # TODO

# minhash/ SBT related
cli.add_command(minhash)    # TODO
cli.add_command(sbt_index)  # TODO

# io
cli.add_command(io)         # TODO

# alignment
cli.add_command(msa_encode)  # TODO
cli.add_command(msa_decode)  # TODO
cli.add_command(sam_encode)  # TODO
cli.add_command(sam_decode)  # TODO

# tree
cli.add_command(tree_encode)  # TODO
cli.add_command(tree_decode)  # TODO
