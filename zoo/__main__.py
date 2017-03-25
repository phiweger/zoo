'''
zoo command line.
'''

import click
import json
from pymongo import MongoClient
from pymongo.errors import DuplicateKeyError


@click.group()
def cli():
    pass


@click.command()
@click.option('--client', default='localhost:27017')
@click.option('--db', default='test')
@click.option('--collection', default='test')
@click.argument('infile', type=click.File('r+'))
def load(infile, client, db, collection):
    '''Load a data cell.'''
    click.echo('Loading data cell.')
    c = MongoClient(client)[db][collection]
    for line in infile:
        try:
            c.insert_one(json.loads(line.strip()))
        except DuplicateKeyError:
            print('Duplicate key, loading aborted.')
            return
    click.echo('# of entries written:', str(c.count()))
    print(c.count())


@click.command()
@click.option('--client', default='localhost:27017')
@click.option('--db', default='test')
@click.option('--collection', default='test')
@click.argument('outfile', type=click.File('w+'))
def dump(outfile, client, db, collection):
    '''Dump a (mongodb) cursor to a data cell.

    For each document, start a new line in the output.

    {"_id":"86853586-5e9...
    {"_id":"689e59b8-514...
    {"_id":"6d9bff35-aab...

    This is important bc/ it circumvents the need to hold more than one record
    in memory, both on import and export. Note also that this is the same
    output format as ...

    $ mongoexport --db foo --collection bar --out bar.json
    ... and can be reimported by
    $ mongoimport --db foo --collection bar2 bar.json
    '''
    click.echo('Dumping data cell.')
    c = MongoClient(client)[db][collection]
    for i in c.find():
        outfile.write(json.dumps(i) + '\n')  # same as zoo.io.dump_json


cli.add_command(load)
cli.add_command(dump)
