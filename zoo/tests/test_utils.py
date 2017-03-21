from copy import deepcopy
import pytest
from zoo.utils import deep_get, deep_set


original = {'a': {'b': 5, 'c': {}, 'd': None, 'e': False, 'f': '', 'g': []}}
d = deepcopy(original)


def test_deep_get():
    assert deep_get(d, 'a.b') == 5


def test_deep_set():
    # Replacing existing values does not work.
    with pytest.raises(TypeError):
        deep_set(d, 'a.b', [1, 2])
        # Key exists, item assignment not allowed w/ replace=False
    deep_set(d, 'a.b', [1, 2], replace=True)

    # Except when they evaluate to False, i.e. [], {}, "", False, None
    deep_set(d, 'a.c', 42)
    deep_set(d, 'a.d', 42)
    deep_set(d, 'a.e', 42)
    deep_set(d, 'a.f', 42)
    deep_set(d, 'a.g', {'foo': 'bar'})

    # We can access and modify the objects in place with deep_get()
    deep_get(d, 'a.b').append(3)
    assert deep_get(d, 'a.b')[2] == 3

    deep_get(d, 'a.g').update({'bar': 'foo'})
    assert deep_get(d, 'a.g')['bar'] == 'foo'

    # We can create new, nested keys.
    with pytest.raises(KeyError):
        deep_set(d, 'a.new.nested.path', 5)
        # 'Key not present. Use "force=True" to create key.'
    deep_set(d, 'a.new.nested.path', 5, force=True)
    assert deep_get(d, 'a.new') == {'nested': {'path': 5}}

    # A new (nested) key can only be created if the "root" is a dict.
    with pytest.raises(AttributeError):
        deep_set(d, 'a.new.nested.path.below', 5, force=True)
        # 'int' object has no attribute 'setdefault'
    deep_set(d, 'a.new.nested.path', {}, replace=True)
    deep_set(d, 'a.new.nested.path.below', 5, force=True)
    assert deep_get(d, 'a.new.nested.path.below') == 5
