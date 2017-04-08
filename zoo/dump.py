import json


'''
TODO: is this quicker/ better?

from bson.json_util import dumps  # stackoverflow, 30333299
with open(path, 'w+') as outjson:
    outjson.write(dumps(db.collection.find(), indent=4))
# fasta is 350 kb, JSON 359
'''


def dump_json(path, cursor):
    '''Export a (mongodb) cursor to JSON format.

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
    with open(path, 'w+') as outjson:
        for i in cursor:
            outjson.write(json.dumps(i) + '\n')


def dump_fasta(query):
    print(next(query))

