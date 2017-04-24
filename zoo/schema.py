import json
from zoo import get_schema


a = {'a': 1, 'b': 2}
b = {'c': 5, 'd': 5}


def schema_combine(schemas, additional=False, required=['_id', 'seq']):
    '''Combine multiple schemas predefined in zoo.

    Example:

    schemas = ['core', 'metadata', 'derivative']
    '''
    d = {}
    for i in 'core metadata relative'.split(' '):
        with open(get_schema(i + '.json'), 'r') as schema:
            d.update(json.load(schema))
    return d




