'''
Notes:

- https://github.com/viehwegerlib/zoo/issues/136
- the "required" field in a JSON schema can be nested, stackoverflow,
  31664705, see also https://github.com/Julian/jsonschema/issues/331
- what about additionalProperties?
- the first level needs to be core, and it cannot coexist in a flattened
  way, i.e. (a,b(c,d)) is not valid (because we'd get key collisions:
  which "description" to use, for example); another way to think about it:
  there can only be one root
- if some schema is nested, it must contain a "properties" field, e.g.
  (a(b,c(d,e))) a and c need to contain the "properties" field
- https://github.com/geraintluff/tv4 in JavaScript
'''


import glob
import json
import re


def find_schema(name, fp):
    '''Given a name and a file path, try to find the schema and import.

    Returns "None" if <name> not found in <fp>.

    This works recursively, i.e. only a top directory needs to be
    provided. If a schema (fragment) is called "name", "name.json" will be
    imported.
    '''
    # only works for Python >= 3.5, stackoverflow, 2186525

    # handle fp "/", i.e. both "path" and "path/" are acceptable
    if fp[-1] != '/':
        fp += '/'

    # note that schema names must be unique for this to work
    for file in glob.iglob(fp + '**/*.json', recursive=True):
        if '/' + name + '.json' in file:
            with open(file, 'r') as schema:
                return(json.load(schema))

    # If function does not return, i.e. name not found.
    raise Exception(
        'Schema (fragment) "' + name + '" not found.\nAbort!')


def combine_schema(tree, fp):
    '''
    Algo (pseudocode):

    init queue
    if "(" append {} to queue
    if name append to []
    if ")" pop last element from queue and insert into penultimate[-1]
    if len(queue) == 1, return its content (i.e. the last remaining element)

    Example:
    import json
    from jsonschema import validate
    from zoo import get_schema
    from zoo.schema import combine_schema

    schema = combine_schema('(a(b,c))', fp=get_schema('test_schemas'))
    print(json.dumps(schema, indent=4, sort_keys=True))
    rec = {'_id': '1', 'zoo': 'le', 'animal': 'lion', 'c':{'very': True}}
    validate(rec, schema)  # fails, bc/ we "danger" is required

    # See what happens when the root is not unique:
    # (a,b,c)
    # ... and when you forget to close parentheses in the tree syntax:
    # (a(b,c)
    # ... and when things are nested that do not contain a "properties" field:
    # (a(b(b,c),c))
    # ... and when schema fragments are specified that are not found:
    # (a(c,d))
    '''
    i = 0
    queue = []  # LIFO queue of dicts

    while i < len(tree):
        # print('i:', i)
        # print('string:', tree[i])
        # print('queue:', queue)

        # 1. case: open backet
        if tree[i] == '(':
            queue.append([])  # {'properties': {}}
            i += 1  # advance string by one character

        # 2. case: closed bracket
        elif tree[i] == ')':
            last = queue.pop()  # like [{'a': {}}, {'b': {}}]
            # print('last', last)
            # print('queue', queue)
            try:
                penultimate = queue[-1]
                # print('penultimate', penultimate)
            except IndexError:  # queue is empty
                if len(last) > 1:  # unpack queue list
                    print('The tree string does not contain a unique root.')
                    print('Abort!')
                    return
                else:
                    result = last[0]
                    key = list(result.keys()).pop()
                    return result[key]
                    # This is where the function returns if successfull.
            # if queue not empty continues here
            key = list(penultimate[-1].keys()).pop()
            for j in last:
                try:
                    penultimate[-1][key]['properties'].update(j)
                except KeyError:
                    print(
                        '''Schemas that nest, e.g. c in c(d,e), must contain a "properties" field.''')
                    print('Abort!')
                    return
            i += 1

        # 3. case: comma
        elif tree[i] == ',':
            i += 1  # advance in tree string

        # 4. case: "leaf label", if we want to use tree terminology
        else:
            leaf = re.match(
                '(\(|,)(\'|\")*([A-Za-z\_\-\|\.0-9\?\/]+)(\'|\"|)(\[)*',
                tree[i-1:i+200]
                )
            if leaf is not None:
                queue[-1].append(
                    {leaf.group(3): find_schema(leaf.group(3), fp)}
                    )
                # TODO: catch key error of case that file not found:
                # it will then simply be containing {}
                i += len(leaf.group(3)) + \
                    leaf.group().count("'") + leaf.group().count('"')

    print('Make sure to provide the right tree syntax.')
    print('For example, did you close all parentheses?')
    print('Abort!')
    return


# DEPRECATED, kept for reference


# def iter_schemas(fps):
#     '''Given an iterator of file paths (JSON), return a generator of hash maps.

#     Example:
#     from zoo.schema import iter_schema
#     names = ['core', 'derivative']
#     fps = [get_schema(i + '.json') for i in names]
#     gen = iter_schemas(fps)

#     '''
#     l = []
#     for i in fps:
#         with open(i, 'r') as schema:
#             l.append(json.load(schema))
#     yield from l


# def combine(schema_iterator, additional=False, required=['_id', 'seq']):
#     '''Combine multiple schemas predefined in zoo.

#     fps         A list of file paths to schema components.

#     required    A list of keys that must be present. If nothing is required,
#                 set to anything "falsy", i.e. that evaluates to False, e.g.
#                 None, [], etc.
#     additional  Are additional properties allowed? Sets the
#                 "additionalProperties" field in the JSON schema.

#     Returns a (JSON) schema, that is like a document template, against which
#     we can then validate documents.

#     Example:

#     from jsonschema import validate
#     from zoo.schema import iter_schemas, combine

#     names = ['core', 'derivative']
#     fps = [get_schema(i + '.json') for i in names]
#     instance = {'_id': '1', 'seq': 'ACTG'}

#     schema = combine(iter_schemas(fps), additional=False, required=None)
#     validate(instance, template)  # returns None, all's well

#     # add a property, which is not one of the allowed properties
#     instance['foo'] = 'bar'
#     validate(instance, schema)  # throws ValidationError
#     '''
#     d = {}

#     if required:
#         d['required'] = required

#     if not additional:
#         d['additionalProperties'] = False

#     for i in schema_iterator:
#         d = update(d, i)

#     # We loose the description, title etc. These fields are merely informative.
#     del d['description']
#     del d['title']
#     return d


# def is_valid(instance, ):
#     '''Validate a hash map against a zoo JSON schemas.

#     zoo compose --fp /optional/path core,derivative(allow,recursion?) template.json
#     # if it does not find them at top directory will recursively enter subfolders, like node_modules
#     zoo init --validate path/to/template.json
#     '''

#     fragments = ['core', 'derivative']
#     schema = combine(fps(fragments), additional=False, required=None)
