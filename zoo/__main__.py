'''
zoo command line.
'''

import click
import json


@click.group()
def cli():
    pass


@click.command()
# @click.option('--count', default=1, help='Number of greetings.')
@click.argument('input', type=click.File('r+'))
def load(input):
    click.echo('Loading the data cell.')
    click.echo(json.load(input))


@click.command()
def dropdb():
    click.echo('Dropped the database')


cli.add_command(load)
cli.add_command(dropdb)
