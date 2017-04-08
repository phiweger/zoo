import json
from pymongo import MongoClient


def env_switch(client, client_env, db, db_env, cell, cell_env):
    # client
    if not client_env:
        client = MongoClient(client)
    else:
        client = MongoClient(client_env)
    # db
    if not db_env:
        database = client[db]
    else:
        database = client[db_env]
    # cell
    if not cell_env:
        collection = database[cell]
    else:
        collection = database[cell_env]
    return collection


def eval_query(collection, query, selection):
    '''Given filepath to JSON, evaluate query and selection, return cursor.'''
    if selection:
        selection = {k: 1 for k in selection.split(',')}
    if query:
        with open(query, 'r+') as file:
            print('Evaluating query.')
            return collection.find(json.load(file), selection)
            # passing None to selection defaults to "return all documents"
    else:
        print('All entries used.')
        return collection.find()  # whole collection


def eval_pipeline(collection, pipeline):
    '''Given filepath to JSON, evaluate pipeline, return cursor.

    Example:

    # pipeline.json
    [
        {
            "$match": {
                "meta.geo.cou": "sierra_leone"
            }
        },
        {
            "$sample": {
                "size": 15
            }
        }
    ]

    '''
    # client = MongoClient('localhost:27017')
    # c = client['ebola']['makona']

    with open('pipeline.json', 'r+') as file:
        return collection.aggregate(json.load(file))


