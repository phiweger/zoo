import click
import json
from progressbar import ProgressBar, UnknownLength
from pymongo import MongoClient
from pymongo.errors import DuplicateKeyError
from sourmash_lib import Estimators, signature, signature_json
from zoo.hash import hash_dict


@click.command()
@click.option('--client', default='localhost:27017')
@click.option('--db')
@click.option('--cell')
@click.argument('file', type=click.File('r+'))
def add(file, client, db, cell):
    '''Load a data cell.

    Example:

    zoo --client localhost:27017 --db zika --cell survey survey.json
    zoo --client localhost:27017 --db zika --cell survey *
    '''
    click.echo('Loading data cell.')
    c = MongoClient(client)[db][cell]
    for line in file:
        try:
            c.insert_one(json.loads(line.strip()))
        except DuplicateKeyError:
            print('Duplicate key, loading aborted.')
            return
    print(c.count(), 'documents inserted in collection', cell + '.\nDone.')


@click.command()
@click.option('--client', default='localhost:27017')
@click.option('--db', required=True)
@click.option(
    '--cell', required=True,
    help='Filename w/o extension.')
@click.option(
    '--ksize', default='16,31',
    help='Comma separated list of k-mer sizes. Larger is more specific.')
@click.option('--n', default=1000)
@click.argument('file', type=click.Path())
def commit(file, client, db, cell, ksize, n):
    '''Dump a (mongodb) cursor to a data cell.

    For each document, start a new line in the output.

    \b
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
    db = MongoClient(client)[db]

    # initialize minhash
    ksize = [int(i) for i in ksize.split(',')]
    dk = {k: Estimators(ksize=k, n=n) for k in ksize}

    bar = ProgressBar(max_value=UnknownLength)
    counter = 0
    with open(file + '.json', 'w+') as f:
        for d in db[cell].find():
            counter += 1

            # calculate fresh md5 hash for each record
            _id = d.pop('_id')
            # Neither the primary key (because it is random)
            # nor the checksum should figure in the checksum.
            try:
                del d['md5']
            except KeyError:
                pass
            d['md5'] = hash_dict(d)
            d['_id'] = _id
            f.write(json.dumps(d, indent=None, sort_keys=True) + '\n')

            # update aggregate minhash for collection
            for v in dk.values():
                v.add_sequence(d['sequence'], force=True)

            # update progress bar
            bar.update(counter)

    # save minhash
    for k, v in dk.items():
        dk.update({
            k: signature.SourmashSignature(
                estimator=v, name=cell, email='', filename='')})
    # print('\n', ksize[0], ksize[1], n)

    with open(file + '.zoo', 'w+') as f:
        signature_json.save_signatures_json(
            dk.values(), fp=f, indent=4, sort_keys=True)
    click.echo('\nDone.')


@click.command()
def diff(outfile, client, db, collection):
    print('Trying.')


@click.command()
def pull():
    # assume there is md5 from commit
    print('Trying.')


@click.command()
def drop():  # cell
    print('Trying.')


@click.command()
def destroy():  # drop database entirely
    print('Trying.')















