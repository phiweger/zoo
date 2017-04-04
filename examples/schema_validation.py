import json
from jsonschema import validate, RefResolver
import os
import re
from zoo import get_schema


# valid
a = 'A/mallard/Interior Alaska/6MP0155/2006(H3N8)'
b = 'A/blue-winged Teal/Minnesota/AI09-2977/2009(H4N8)'
c = 'A/Peru/PER175/2012(H3N2)'

# not valid

# missing components
d = 'A/mallard/Interior Alaska/6MP0155/2006'
# year invalid
e = 'A/Peru/PER175/4912(H3N2)'
# appended information
f = 'A/Peru/PER175/2012(H3N2)-unsure'
# appended information
g = 'A/Peru/PER175/2012(HXN2)'

# stackoverflow, 4374185
p = re.compile('^A/[.+?/]?.+?/.+?/(19|20)\d{2}\(H\dN\d\)$')
# in JSON format, the pattern needs a slight change: escaping needs 2 "\"
# '^A/[.+?/]?.+?/.+?/(19|20)\\d{2}\\(H\\dN\\d\\)$'


for i in 'abcdefg':
    try:
        print('match:', re.match(p, eval(i)).group(0))
    except AttributeError:  # 'NoneType' object has no attribute 'group'
        print('invalid string:', eval(i))
'''
match: A/mallard/Interior Alaska/6MP0155/2006(H3N8)
match: A/blue-winged Teal/Minnesota/AI09-2977/2009(H4N8)
match: A/Peru/PER175/2012(H3N2)
invalid string: A/mallard/Interior Alaska/6MP0155/2006
invalid string: A/Peru/PER175/4912(H3N2)
invalid string: A/Peru/PER175/2012(H3N2)-unsure
invalid string: A/Peru/PER175/2012(HXN2)
'''


'''
Schema validation.
'''


# mock schema
with open(get_schema('specific/influenza_a.json'), 'r+') as file:
    schema = json.load(file)


valid = 'A/Peru/PER175/2012(H3N2)'
invalid = 'A/Peru/PER175/9999(H3N2)'

i = {
    '_id': '1',
    'sequence': 'ATG...',
    'segment': None,
    'nomenclature': valid,  # testing pattern matching
    'annotation': {'id': 4}  # testing "$ref" JSON pointer to another fragment
    }

j = {
    '_id': '1',
    'sequence': 'ATG...',
    'segment': None,
    'nomenclature': invalid
    }


'''
For some reason JSON pointer to JSON fragment in directory invalid, like:
{"$ref": "fragments/test_extension.json#"}

Known issue: https://github.com/Julian/jsonschema/issues/313

However, this syntax is allowed, see:
https://spacetelescope.github.io/understanding-json-schema/structuring.html#structuring
'''

# workaround, #313
sSchemaDir = os.path.dirname(os.path.abspath(get_schema('core.json')))
oResolver = RefResolver(
    base_uri='file://' + sSchemaDir + '/', referrer=schema)

validate(i, schema=schema, resolver=oResolver)
# returns None, i.e. success

validate(j, schema=schema, resolver=oResolver)
'''
ValidationError: 'A/Peru/PER175/9999(H3N2)' does not match
'^A/[.+?/]?.+?/.+?/(19|20)\\d{2}\\(H\\dN\\d\\)$'

Failed validating 'pattern' in schema['properties']['nomenclature']:
    {'pattern': '^A/[.+?/]?.+?/.+?/(19|20)\\d{2}\\(H\\dN\\d\\)$',
     'type': 'string'}

On instance['nomenclature']:
    'A/Peru/PER175/9999(H3N2)'
'''
