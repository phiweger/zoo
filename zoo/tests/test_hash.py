import json
from zoo.utils import ordered
from zoo.hash import hash_content


a = json.loads("""
{
    "errors": [
        {"error": "invalid", "field": "email"},
        {"field": null, "error": "required", "lala": {"new": 4}}
    ],
    "list": [5,3,2,4,1,6],
    "success": false
}
""")

b = json.loads("""
{
    "success": false,
    "errors": [
        {"lala": {"new": 4}, "error": "required", "field": null},
        {"error": "invalid", "field": "email"}
    ],
    "list": [1,2,3,4,5,6]
}
""")


def test_ordered():
    assert ordered(a) == ordered(b)


def test_content_hash():
    assert hash_content(a, 'md5') == hash_content(b, 'md5')
