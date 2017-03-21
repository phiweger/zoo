from copy import deepcopy
from zoo.utils import deep_get, deep_set


original = {'a': {'b': 5}}
d = deepcopy(original)


def test_deep_get():
    assert deep_get(d, 'a.b') == 5

    deep_set(d, 'a.b', [1, 2])
    deep_get(d, 'a.b').append(3)

    assert deep_get(d, 'a.b')[2] == 3



def test_deep_set():
    # modify int to {} and then force create nested key there

# make sure original has not been modified
# assert original == {'a': {'b': 5}}