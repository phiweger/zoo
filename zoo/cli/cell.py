import click
import json
from progressbar import ProgressBar, UnknownLength
from pymongo import MongoClient
from pymongo.errors import DuplicateKeyError
from sourmash_lib import Estimators, signature, signature_json
from uuid import uuid4
from zoo.hash import hash_dict
from zoo.utils import deep_get


'''
Helper functions.
'''


def abort_if_false(ctx, param, value):
    '''http://click.pocoo.org/5/options/#prompting'''
    if not value:
        ctx.abort()


'''
CLI
'''


@click.command()
@click.option('--client', default='localhost:27017')
@click.option('--db', required=True)
@click.option('--cell', required=True, help='Cell name.')
@click.option(
    '--primkey', default='_id',
    help='What is the primary key to judge duplicates by? UUID, GenBank?')
@click.argument('file', type=click.File('r+'))
def add(file, client, db, cell, primkey):
    '''Load a data cell.

    An alternative primary key can be specified to insert documents. This
    is useful in the case where the data cell comes from a collaborator
    who uses a different set of UUIDs as we do. In this case, these identifiers
    do not reflect, whether an entry is a duplicate.

    Example:

    \b
    $ zoo add --client localhost:27017 --db zika --cell t5 zoo/data/cell_a.json
    Loading data cell.
    3 documents inserted in collection t5.
    0 duplicates skipped.
    Done.

    \b
    $ zoo add --db zika --cell t5 --primkey genbank.a zoo/data/cell_b.json
    Loading data cell.
    Index created on field "genbank.a".
    1 documents inserted in collection t5.
    3 duplicates skipped.
    Done.
    '''
    click.echo('Loading data cell.')
    c = MongoClient(client)[db][cell]
    duplicates = 0
    if primkey == '_id':
        for line in file:
            try:
                c.insert_one(json.loads(line.strip()))
            except DuplicateKeyError:
                duplicates += 1
                pass
    else:
        # index primkey if it does not exists yet
        if primkey not in c.index_information():
            c.create_index(primkey, unique=True, name=primkey)
            print('Index created on field', '"' + primkey + '".')
        for line in file:
            d = json.loads(line.strip())
            if c.find_one({primkey: deep_get(d, primkey)}):  # no duplicate
                duplicates += 1
            else:
                c.insert_one(d)

    print(
        c.count() - duplicates, 'documents inserted in cell', cell + '.')
    if duplicates > 0:
        print(duplicates, 'duplicates skipped.\nDone.')


@click.command()
@click.option('--client', default='localhost:27017')
@click.option('--db', required=True)
@click.option(
    '--cell', required=True,
    help='Cell name.')
@click.option(
    '--ksize', default='16,31',
    help='Comma separated list of k-mer sizes. Larger is more specific.')
@click.option('--n', default=1000)
@click.argument(
    'file', type=click.Path())
def commit(file, client, db, cell, ksize, n):
    '''Dump a (mongodb) cursor to a data cell.

    For each document, start a new line in the output.

    file argument: Filename prefix w/o extension.

    \b
    {"_id":"86853586-5e9...
    {"_id":"689e59b8-514...
    {"_id":"6d9bff35-aab...

    This is important bc/ it circumvents the need to hold more than one record
    in memory, both on import and export. Note also that this is the same
    output format as ...

    \b
    $ mongoexport --db foo --collection bar --out bar.json
    ... and can be reimported by
    $ mongoimport --db foo --collection bar2 bar.json

    Example:

    $ zoo commit --db zika --cell survey --n 5 surveytest
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
    '''
    for each entry in pulled file:
    1. calculate fresh md5
    2. compare to existing
    3. if != drop entry end replace with pulled one
    4. we can assume primary keys same, because we pull changes from existing.
    '''
    print('Trying.')


@click.command()
def status():
    '''
    \b
    - list cells, num_entries
    - verbose: find_one() in each cell but truncate sequence field before print
    - include .zoo metadata in the report in the future
    '''
    print('Trying.')


@click.option('--client', default='localhost:27017')
@click.option('--db', required=True)
@click.option('--cell', required=True, help='Cell name.')
@click.option('--force', is_flag=True, callback=abort_if_false,
              expose_value=False,
              prompt='Are you sure you want to drop the db?')
@click.command()
def drop(client, db, cell):  # cell
    '''Delete a cell.

    If a cell did not exist, it is dropped nevertheless. Ask your favourite
    philosopher what it means to delete something that did not exist.

    Example:

    \b
    $ zoo drop --db zika --cell survey
    Are you sure you want to drop the db? [y/N]: y
    Dropped cell "animals" from database "zika".

    \b
    $ zoo drop --db zika --cell survey --force  # no confirmation
    Dropped cell "animals" from database "zika".
    '''
    c = MongoClient(client)[db][cell]
    c.drop()
    print('Dropped cell', '"' + cell + '"', 'from database', '"' + db + '".')


@click.command()
def destroy():  # drop database entirely
    print('Trying.')


@click.option('--client', default='localhost:27017')
@click.option('--db', required=True)
@click.option('--cell', required=True, help='Cell name.')
@click.argument('file', type=click.Path())
@click.command()
def init(file, client, db, cell):  # load json to mongodb and assign UUID
    '''
    If we write, no file extension needed, if we read, needed to indicate file.
    See also "zoo add ..." (read) vs. "zoo commit" (write).

    Example:

    \b
    $ zoo init --db zika --cell animals zoo/data/cell_a.json
    Initializing data cell.
    Inserted 3 entries into "animals".
    '''
    click.echo('Initializing data cell.')
    c = MongoClient(client)[db][cell]
    with open(file, 'r+') as f:
        for line in f:
            d = json.loads(line.strip())
            d['_id'] = str(uuid4())
            c.insert_one(d)
    print('Inserted', c.count(), 'entries into cell', '"' + cell + '".')














