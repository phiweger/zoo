'''
Utils for Machine Learning.
'''


def sample(collection, field, stratify=True, n):
    '''

    '''
    pass


client = MongoClient("localhost:27017")
db = client["zoo"]
col = db.get_collection('influenza_a_virus')


pipeline = [ 
    { '$group': { '_id': "metadata.host"}  },
    { '$group': { '_id': 1, 'count': { '$sum': 1 } } }
];

//
// Run the aggregation command
//
R = db.command( 
    {
    "aggregate": "influenza_a_virus" , 
    "pipeline": pipeline
    }
);
printjson(R);



q = col.find()

birds = col.aggregate([
{'$match': {
    'annotation.name': 'HA',
    'metadata.host': 'Avian'
    }},
{'$sample': {
    'size': 15
    }}
])


col.distinct('metadata.host')






