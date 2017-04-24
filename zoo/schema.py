import json
from zoo import get_schema
from zoo.utils import update


def fps(names):
    '''Get the filepaths of predefined schemas in zoo.

    Example:
    from zoo.schema import fps
    fps('core metadata relative'.split(' '))

    '''
    return [get_schema(i + '.json') for i in names]


def combine(fps, additional=False, required=['_id', 'seq']):
    '''Combine multiple schemas predefined in zoo.

    fps         A list of file paths to schema components.

    required    A list of keys that must be present. If nothing is required,
                set to anything "falsy", i.e. that evaluates to False, e.g.
                None, [], etc.
    additional  Are additional properties allowed? Sets the
                "additionalProperties" field in the JSON schema.

    Returns a (JSON) schema, that is like a document template, against which
    we can then validate documents.

    Example:

    from jsonschema import validate
    from zoo.schema import fps, combine

    instance = {'_id': '1', 'seq': 'ACTG'}
    fragments = ['core', 'derivative']
    schema = combine(fps(fragments), additional=False, required=None)
    validate(instance, template)  # returns None, all's well

    # add a property, which is not one of the allowed properties
    instance['foo'] = 'bar'
    validate(instance, template)  # throws ValidationError
    '''
    d = {}
    for i in fps:
        with open(i, 'r') as schema:
            d = update(d, json.load(schema))
            # print(d)

    if required:
        d['required'] = required

    if not additional:
        d['additionalProperties'] = False
    return d
