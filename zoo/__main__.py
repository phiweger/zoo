'''
zoo command line.
'''

import click
from .cli.cell import init, add, commit, diff, pull, drop, destroy


@click.group()
def cli():
    pass


cli.add_command(init)
cli.add_command(add)
cli.add_command(commit)
cli.add_command(diff)
cli.add_command(pull)
cli.add_command(drop)
cli.add_command(destroy)
