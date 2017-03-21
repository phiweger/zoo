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


'''
In [57]: d = {'a': {'b': ''}}

In [58]: deep_set(d, 'a.b', 5)

In [59]: d
Out[59]: {'a': {'b': 5}}

In [60]: deep_set(d, 'a.b', '')
---------------------------------------------------------------------------
TypeError                                 Traceback (most recent call last)
<ipython-input-60-7f2b9aed0fea> in <module>()
----> 1 deep_set(d, 'a.b', '')

<ipython-input-56-2d2df156a8bd> in deep_set(d, key, value, force, replace)
     51             # i.e. field is non-False and we don't want to replace it
     52             raise TypeError(
---> 53                 'Key exists, item assignment not allowed w/ replace=False')
     54         else:
     55             reduce(dict.__getitem__, keys, d)[latest] = value

TypeError: Key exists, item assignment not allowed w/ replace=False

In [61]: deep_set(d, 'a.b', '', replace=True)

In [62]: d
Out[62]: {'a': {'b': ''}}

In [63]: deep_set(d, 'a.b.c', 5)
---------------------------------------------------------------------------
KeyError                                  Traceback (most recent call last)
<ipython-input-63-0c6930b96253> in <module>()
----> 1 deep_set(d, 'a.b.c', 5)

<ipython-input-56-2d2df156a8bd> in deep_set(d, key, value, force, replace)
     60         if not force:
     61             raise KeyError(
---> 62                 'Key not present. Use "force=True" to create key.')
     63         else:
     64             for k in keys:

KeyError: 'Key not present. Use "force=True" to create key.'

In [64]: deep_set(d, 'a.b.c', 5, force=True)
---------------------------------------------------------------------------
AttributeError                            Traceback (most recent call last)
<ipython-input-64-185110f2192f> in <module>()
----> 1 deep_set(d, 'a.b.c', 5, force=True)

<ipython-input-56-2d2df156a8bd> in deep_set(d, key, value, force, replace)
     64             for k in keys:
     65                 d = d.setdefault(k, {})
---> 66             d.setdefault(latest, value)
     67

AttributeError: 'str' object has no attribute 'setdefault'

In [65]: d
Out[65]: {'a': {'b': ''}}

In [66]: deep_set(d, 'a.b', {})

In [67]: deep_set(d, 'a.b.c', 5, force=True)

In [68]: d
Out[68]: {'a': {'b': {'c': 5}}}
'''