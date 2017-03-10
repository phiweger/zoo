import numpy as np


def entropy_select(arr, n):
    '''
    # Do negative control by shuffling arr

    {0, ..., 1}: top proportion
    >1: top number of position

    Usage:
    a = np.array([3, 4, 8, 23, 1, 43])
    entropy_select(a, 3)
    # Will select top 3 elements
    # array([2, 3, 5])
    '''
    if n < 1:
        print('Will select top', int(n * 100), 'percent.')
        threshold = np.percentile(arr, (1 - n) * 100)
        return np.where(arr > threshold)[0]
    elif n >= 1:
        print('Will select top', int(n), 'elements.')
        return arr.argsort()[-int(n):]  # stackoverflow, 6910641


# '''
# Utils for Machine Learning.
# '''

# from sklearn.model_selection import train_test_split
# import pymongo


# # match, group, sample

# '''
# or match group then match sample (because group loses doc info?)


# def(match=d1, group=d2, sample=15 or "stratified" or weight list/dict)
# '''


# hosts = 'Human Avian Swine'.split(' ')
# match = {
#     'annotation.name': 'HA',
#     'metadata.host': {'$in': hosts},
#     'metadata.date.y': {'$gte': 2000}
#     }


# q = collection.find(match)
# gen = ((i['_id'], i['metadata']['host']) for i in q)
# df = pd.DataFrame.from_records(gen)




# group = {
#     '_id': '$metadata.host',
#     'cnt': {'$sum': 1}
#     }

# # group = {'$group': {
# #     '_id': {
# #         'host': '$metadata.host',
# #         'country': '$metadata.country'
# #         },
# #     'cnt': {'$sum': 1}
# #     }}

# sample = {'$sample': {
#     'size': n
#     }}

# # pipe1 aggregates "what's there", so select number to sample
# pipe1 = [{'$match': match}, {'$group': group}]
# q = col.aggregate(pipe1)
# d = {}
# for i in q:
#     d[i['_id']] = i['cnt']
# # {'Avian': 9960, 'Human': 14811, 'Swine': 3136}


# gen = ((i['_id'], i['metadata.host']) for i in q)



# selection = {'Avian': 10000, 'Human': 15000, 'Swine': 3359}
# weights = {'test': 0.3, 'train': 0.5, 'validate': 0.2}  # assert sum to 1








# def test_train_split():
#     '''
#     *arrays : sequence of indexables with same length / shape[0]

#     Allowed inputs are lists, numpy arrays, scipy-sparse matrices or pandas dataframes.

# test_size : float, int, or None (default is None)

#     If float, should be between 0.0 and 1.0 and represent the proportion of the dataset to include in the test split. If int, represents the absolute number of test samples. If None, the value is automatically set to the complement of the train size. If train size is also None, test size is set to 0.25.

# train_size : float, int, or None (default is None)

#     If float, should be between 0.0 and 1.0 and represent the proportion of the dataset to include in the train split. If int, represents the absolute number of train samples. If None, the value is automatically set to the complement of the test size.

# random_state : int or RandomState

#     Pseudo-random number generator state used for random sampling.

# stratify : array-like or None (default is None)

#     If not None, data is split in a stratified fashion, using this as the class labels.

#     '''


# match = {'$match': {
#     'annotation.name': 'HA',
#     'metadata.host': 'Avian'
# }}

# sample = {'$sample': {
#     'size': 30
# }}

# pipe2 = [match, sample]
# q = col.aggregate(pipe2)



# def bucket(cursor, test, seed=None):
#     if seed:
#         random.seed(seed)
#     r = random.uniform(0, 1)
#     for key in weights.keys():



# '''
# (id, label) array/ dataframe
# test_train_split():
# '''



# split(collection, match, selection, weights):
#     for tag, n in selection.items():
    
#         match = {'$match': {
#             'annotation.name': 'HA',
#             'metadata.host': tag
#         }}
    
#         sample = {'$sample': {
#             'size': n
#         }}
    
#         pipe2 = [match, sample]
#         q = col.aggregate(pipe2)

# '''
# The sampling in mongoDB is without replacement. However, if a write is
# performed on the database while sampling, a (modified) document may be
# returned more than once. See "cursor isolation" in docs:

# > As a cursor returns documents, other operations may interleave with the
# query. For the MMAPv1 storage engine, intervening write operations on a
# document may result in a cursor that returns a document more than once
# if that document has changed. To handle this situation, see the information
# on snapshot mode.
# '''



# def group_one(collection, field):
#     '''
#     syntactic sugar

#     group(col, 'metadata.host')
#     group(col, 'metadata.date.y')

#     {'': 224,
#      'Avian': 91698,
#      'Bat': 32,
#      'Camel': 8,
#      ...}
#     '''
#     pipeline = [{'$group': {'_id': '$' + field, 'cnt': {'$sum': 1}}}]
#     q = collection.aggregate(pipeline)
#     d = {}
#     for i in q:
#         d[i['_id']] = i['cnt']
#     return d


# def group_many(collection, field):
#     '''
#     syntactic sugar

#     group_many(col, ['metadata.host', 'metadata.country', 'metadata.date.y'])

#     TODO: allow time range grouping: stackoverflow, 27750974
#     '''
#     d = {}
#     for i in field:
#         key = i.split('.')[-1]
#         d[key] = '$' + i

#     # https://gist.github.com/clarkenheim/fa0f9e5400412b6a0f9d
#     pipeline = [{
#         '$group': {
#             '_id': d,
#             'cnt': {'$sum': 1}
#             }}]
#     try:
#         q = collection.aggregate(pipeline)
#     except pymongo.errors.OperationFailure:
#         print('Check if <field> argument is a list.')
#         return 1

#     e = {}
#     for i in q:
#         e[tuple(i['_id'].values())] = i['cnt']
#     return e


# g = group_one(col, 'metadata.host')
# d = {host: g[host] for host in hosts}


# # sample_n
# # sample_ml
# sample(match, group, n):
#     pipeline1 = [xmatch, xgroup]
#     pipeline2 = [xmatch, n]

#     # sample and randomly assign to either:
#     return test, validation, train

# # d
# {test: [{document_1}, {2}, ...], train: [{}, {}, ...], validate: [{}, {}, ...]}


# split(match, group, weights):
#     '''
#     '''
#     pass


# match = {
#     'annotation.name': 'HA',
#     'metadata.host': {'$in': hosts}
#     }

# # group = single or many calling appropriate function

# group = {
#         '_id': {
#             'host': '$metadata.host',
#             'country': '$metadata.country'},
#         'cnt': {'$sum': 1}
#         }

# n = 15


# def sample(q, n)


# client = MongoClient("localhost:27017")
# db = client["zoo"]
# col = db.get_collection('influenza_a_virus')


# n = 15
# hosts = 'Avian Swine Human'.split(' ')


# available(col, 'metadata.host', hosts, )


# # doesn't work on documents
# # so match - group
# # same match - sample
# # stackoverflow, 11782566
# pipeline = [
#     {'$match': {
#         'annotation.name': 'HA',
#         'metadata.host': {'$in': hosts}
#         }},
#     {'$group': {
#         '_id': {
#             'host': '$metadata.host',
#             'country': '$metadata.country'},
#         'cnt': {'$sum': 1}
#         }},
#     {'$sample': {
#         'size': n
#         }}
#         ]

# birds = col.aggregate(pipeline)


# col.distinct('metadata.host')

















