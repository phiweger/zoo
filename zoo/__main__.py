'''
zoo command line.
'''

import click
import ijson


@click.group()
def cli():
    pass


@click.command()
# @click.option('--count', default=1, help='Number of greetings.')
@click.argument('input', type=click.File('r+'))
def load(input):
    click.echo('Loading the data cell.')
    click.echo(json.load(input))
    # db[args.collection].insert(json.loads(line))
    # http://trimc-db.blogspot.de/2015/01/using-python-driver-in-mongodb-for.html

'''
We have to assume JSON input is larger than memory. We'll use "ijson".
- https://www.dataquest.io/blog/python-json-tutorial/
- https://pypi.python.org/pypi/ijson/

filename = "md_traffic.json"
with open(filename, 'r') as f:
    objects = ijson.items(f, 'meta.view.columns.item')
    columns = list(objects)


with open(fn, 'r+') as input:
    for item in islice(ijson.items(input, 'item'), 5):
    # mongodb dump is a list of JSON files
    # stackoverflow, 19996401
        print(item['metadata'])
    

'''


@click.command()
def dropdb():
    click.echo('Dropped the database')


cli.add_command(load)
cli.add_command(dropdb)
