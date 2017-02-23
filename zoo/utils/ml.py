'''
Utils for Machine Learning.
'''


import pymongo


def group_one(collection, field):
    '''
    group(col, 'metadata.host')
    group(col, 'metadata.date.y')

    {'': 224,
     'Avian': 91698,
     'Bat': 32,
     'Camel': 8,
     ...}
    '''
    pipeline = [{'$group': {'_id': '$' + field, 'cnt': {'$sum': 1}}}]
    q = collection.aggregate(pipeline)
    d = {}
    for i in q:
        d[i['_id']] = i['cnt']
    return d


def group_many(collection, field):
    '''
    group_many(col, ['metadata.host', 'metadata.country', 'metadata.date.y'])

    TODO: allow time range grouping: stackoverflow, 27750974
    '''
    d = {}
    for i in field:
        key = i.split('.')[-1]
        d[key] = '$' + i

    # https://gist.github.com/clarkenheim/fa0f9e5400412b6a0f9d
    pipeline = [{
        '$group': {
            '_id': d,
            'cnt': {'$sum': 1}
            }}]
    try:
        q = collection.aggregate(pipeline)
    except pymongo.errors.OperationFailure:
        print('Check if <field> argument is a list.')
        return 1

    e = {}
    for i in q:
        e[tuple(i['_id'].values())] = i['cnt']
    return e











def sample(q, n)


client = MongoClient("localhost:27017")
db = client["zoo"]
col = db.get_collection('influenza_a_virus')


n = 15
hosts = 'Avian Swine Human'.split(' ')


available(col, 'metadata.host', hosts, )



# stackoverflow, 11782566
pipeline = [
    {'$match': {
        'annotation.name': 'HA',
        'metadata.host': {'$in': hosts}
        }},
    {'$sample': {
        'size': n
        }}
        ]

birds = col.aggregate(pipeline)


col.distinct('metadata.host')

















