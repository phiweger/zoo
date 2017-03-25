'''
zoo command line.
'''

import click


@click.group()
def cli():
    pass


@click.command()
@click.option('--count', default=1, help='Number of greetings.')
def initdb():
    click.echo('Initialized the database')
    for i in range(count):
        print('foo')


@click.command()
def dropdb():
    click.echo('Dropped the database')


cli.add_command(initdb)
cli.add_command(dropdb)
