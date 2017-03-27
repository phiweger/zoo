'''
zoo command line.
'''

import click
from .cli_cell import add, commit, diff, pull


@click.group()
def cli():
    pass


cli.add_command(add)
cli.add_command(commit)
cli.add_command(diff)
cli.add_command(pull)

