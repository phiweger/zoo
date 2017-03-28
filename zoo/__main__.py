'''
zoo command line.
'''

import click
from .cli.cell import init, add, commit, diff, pull, status, drop, destroy


@click.group()
def cli():
    pass


cli.add_command(init)
cli.add_command(add)
cli.add_command(commit)
cli.add_command(diff)       # TODO
cli.add_command(pull)       # TODO
cli.add_command(status)     # TODO
cli.add_command(drop)       # TODO
cli.add_command(destroy)    # TODO


