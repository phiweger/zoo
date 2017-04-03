import json
from jsonschema import validate, RefResolver
import os
from zoo import get_schema


with open(get_schema('core.json'), 'r+') as file:
    schema = json.load(file)


i = {'_id': '1', 'sequence': 'AACTA', 'annotation': {'id': 1}}

'''

for some reason JSON pointer invalid, see
https://github.com/Julian/jsonschema/issues/313

this is deemed acceptible here:
https://spacetelescope.github.io/understanding-json-schema/structuring.html#structuring

{
    "$schema": "http://json-schema.org/draft-04/schema#",
    "title": "Sequence",
    "description": "The core sequence document.",
    "type": "object",
    "properties": {
        "_id" : {"type" : "string"},
        "sequence" : {"type" : "string"},
        "annotation": {"$ref": "testlego.json#"}
    }
}
'''

# workaround, #313
sSchemaDir = os.path.dirname(os.path.abspath(get_schema('core.json')))
oResolver = RefResolver(
    base_uri='file://' + sSchemaDir + '/', referrer=schema)

validate(i, schema=schema, resolver=oResolver)
# returns None, i.e. success
