from jsonschema import validate
from jsonschema.exceptions import ValidationError
import pytest
from zoo import get_schema
from zoo.schema import combine_schema


def test_combine_schema():
    result = combine_schema('(a(b,c))', fp=get_schema('test_schemas'))
    assert result['properties']['c']['required'] == ['danger']

    # Schema fragments are specified but cannot be found.
    with pytest.raises(Exception):
        combine_schema('(a(b,d))', fp=get_schema('test_schemas'))

    # No root specified.
    result = combine_schema('(a,b,c))', fp=get_schema('test_schemas'))
    assert result is None

    # No "properties" field in schema fragment that is intended to nest.
    result = combine_schema('(a,b(a,c))', fp=get_schema('test_schemas'))
    assert result is None

    # Parentheses not closed.
    result = combine_schema('(a(b,c)', fp=get_schema('test_schemas'))
    assert result is None


def test_validation():
    schema = combine_schema('(a(b,c))', fp=get_schema('test_schemas'))
    rec = {'_id': '1', 'zoo': 'le', 'animal': 'lion', 'c': {'very': True}}
    with pytest.raises(ValidationError):
        validate(rec, schema)  # fails, bc/ we "danger" is required

    rec = {'_id': '1', 'zoo': 'le', 'animal': 'lion', 'c': {'danger': True}}
    assert validate(rec, schema) is None
