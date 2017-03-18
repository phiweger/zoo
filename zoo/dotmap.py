'''
This beautiful class has been adopted from:

https://github.com/drgrib/dotmap

... and was slightly modified so as to allow attribute queries of the form:

doc.dget_attr('a.b.c.d')

... w/ a corresponding setter method.
'''


from collections import OrderedDict, MutableMapping
from json import dumps
import functools
from pprint import pprint
from inspect import ismethod


class DotMap(MutableMapping, OrderedDict):
    def __init__(self, *args, **kwargs):
        self._map = OrderedDict()
        self._dynamic = True
        if kwargs:
            if '_dynamic' in kwargs:
                self._dynamic = kwargs['_dynamic']
        if args:
            d = args[0]
            if isinstance(d, dict):
                for k, v in self.__call_items(d):
                    if isinstance(v, dict):
                        v = DotMap(v, _dynamic=self._dynamic)
                    if type(v) is list:
                        l = []
                        for i in v:
                            n = i
                            if type(i) is dict:
                                n = DotMap(i, _dynamic=self._dynamic)
                            l.append(n)
                        v = l
                    self._map[k] = v
        if kwargs:
            for k, v in self.__call_items(kwargs):
                if k is not '_dynamic':
                    self._map[k] = v

    def dget_attr(self, attr, default=object()):  # d .. document
        '''
        modelled after example: stackoverflow, 31174295, much reduced

        Example:

        doc.dget_attr('a.d')  # same as doc.a.d
        # [{'hello': 'A'}, {'hello': 'B'}]
        [i['hello'] for i in doc.dget_attr('a.d')]
        # ['A', 'B']
        '''

        if default is object():
            _getattr = getattr
        else:
            def _getattr(obj, name):
                return getattr(obj, name, default)
        return functools.reduce(getattr, [self] + attr.split('.'))

    def dset_attr(self, attr, val):
        '''
        doc.dget_attr('a.c.d')
        # 4
        doc.dset_attr('a.c.d', 5)
        doc.dget_attr('a.c.d')
        # 5
        '''
        pre, _, post = attr.rpartition('.')
        # 'a.c.d'.rpartition('.') == ('a.c', '.', 'd')
        # .split() does not work for json depth > 2
        return setattr(self.dget_attr(pre) if pre else self, post, val)

    def __call_items(self, obj):
        if hasattr(obj, 'iteritems') and ismethod(getattr(obj, 'iteritems')):
            return obj.iteritems()
        else:
            return obj.items()

    def items(self):
        return self.iteritems()

    def iteritems(self):
        return self.__call_items(self._map)

    def __iter__(self):
        return self._map.__iter__()

    def next(self):
        return self._map.next()

    def __setitem__(self, k, v):
        self._map[k] = v

    def __getitem__(self, k):
        if k not in self._map and self._dynamic and k != '_ipython_canary_method_should_not_exist_':
            # automatically extend to new DotMap
            self[k] = DotMap()
        return self._map[k]

    def __setattr__(self, k, v):
        if k in {'_map', '_dynamic', '_ipython_canary_method_should_not_exist_'}:
            super(DotMap, self).__setattr__(k,v)
        else:
            self[k] = v

    def __getattr__(self, k):
        if k == {'_map', '_dynamic', '_ipython_canary_method_should_not_exist_'}:
            super(DotMap, self).__getattr__(k)
        else:
            return self[k]

    def __delattr__(self, key):
        return self._map.__delitem__(key)

    def __contains__(self, k):
        return self._map.__contains__(k)

    def __str__(self):
        items = []
        for k, v in self.__call_items(self._map):
            # bizarre recursive assignment situation (why someone would do
            # this is beyond me)
            if id(v) == id(self):
                items.append('{0}=DotMap(...)'.format(k))
            else:
                items.append('{0}={1}'.format(k, repr(v)))
        joined = ', '.join(items)
        out = '{0}({1})'.format(self.__class__.__name__, joined)
        return out

    def __repr__(self):
        return str(self)

    def toDict(self):
        d = {}
        for k, v in self.items():
            if type(v) is DotMap:
                # bizarre recursive assignment support
                if id(v) == id(self):
                    v = d
                else:
                    v = v.toDict()
            elif type(v) is list:
                l = []
                for i in v:
                    n = i
                    if type(i) is DotMap:
                        n = i.toDict()
                    l.append(n)
                v = l
            d[k] = v
        return d

    def pprint(self, pformat='dict'):
        if pformat == 'json':
            print(dumps(self.toDict(), indent=4, sort_keys=True))
        else:
            pprint(self.toDict())

    def empty(self):
        return (not any(self))

    # proper dict subclassing
    def values(self):
        return self._map.values()

    # ipython support
    def __dir__(self):
        return self.keys()

    @classmethod
    def parseOther(self, other):
        if type(other) is DotMap:
            return other._map
        else:
            return other

    def __cmp__(self, other):
        other = DotMap.parseOther(other)
        return self._map.__cmp__(other)

    def __eq__(self, other):
        other = DotMap.parseOther(other)
        if not isinstance(other, dict):
            return False
        return self._map.__eq__(other)

    def __ge__(self, other):
        other = DotMap.parseOther(other)
        return self._map.__ge__(other)

    def __gt__(self, other):
        other = DotMap.parseOther(other)
        return self._map.__gt__(other)

    def __le__(self, other):
        other = DotMap.parseOther(other)
        return self._map.__le__(other)

    def __lt__(self, other):
        other = DotMap.parseOther(other)
        return self._map.__lt__(other)

    def __ne__(self, other):
        other = DotMap.parseOther(other)
        return self._map.__ne__(other)

    def __delitem__(self, key):
        return self._map.__delitem__(key)

    def __len__(self):
        return self._map.__len__()

    def clear(self):
        self._map.clear()

    def copy(self):
        return DotMap(self)

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo=None):
        return self.copy()

    def get(self, key, default=None):
        return self._map.get(key, default)

    def has_key(self, key):
        return key in self._map

    def iterkeys(self):
        return self._map.iterkeys()

    def itervalues(self):
        return self._map.itervalues()

    def keys(self):
        return self._map.keys()

    def pop(self, key, default=None):
        return self._map.pop(key, default)

    def popitem(self):
        return self._map.popitem()

    def setdefault(self, key, default=None):
        self._map.setdefault(key, default)

    def update(self, *args, **kwargs):
        if len(args) != 0:
            self._map.update(*args)
        self._map.update(kwargs)

    def viewitems(self):
        return self._map.viewitems()

    def viewkeys(self):
        return self._map.viewkeys()

    def viewvalues(self):
        return self._map.viewvalues()

    @classmethod
    def fromkeys(cls, seq, value=None):
        d = DotMap()
        d._map = OrderedDict.fromkeys(seq, value)
        return d

    def __getstate__(self): return self.__dict__

    def __setstate__(self, d): self.__dict__.update(d)
    # bannerStr

    def _getListStr(self, items):
        out = '['
        mid = ''
        for i in items:
            mid += '  {}\n'.format(i)
        if mid != '':
            mid = '\n' + mid
        out += mid
        out += ']'
        return out

    def _getValueStr(self, k, v):
        outV = v
        multiLine = len(str(v).split('\n')) > 1
        if multiLine:
            # push to next line
            outV = '\n' + v
        if type(v) is list:
            outV = self._getListStr(v)
        out = '{} {}'.format(k, outV)
        return out

    def _getSubMapDotList(self, pre, name, subMap):
        outList = []
        if pre == '':
            pre = name
        else:
            pre = '{}.{}'.format(pre, name)

        def stamp(pre, k, v):
            valStr = self._getValueStr(k, v)
            return '{}.{}'.format(pre, valStr)
        for k, v in subMap.items():
            if isinstance(v, DotMap) and v != DotMap():
                subList = self._getSubMapDotList(pre, k, v)
                outList.extend(subList)
            else:
                outList.append(stamp(pre, k, v))
        return outList

    def _getSubMapStr(self, name, subMap):
        outList = ['== {} =='.format(name)]
        for k, v in subMap.items():
            if isinstance(v, DotMap) and v != DotMap():
                # break down to dots
                subList = self._getSubMapDotList('', k, v)
                # add the divit
                # subList = ['> {}'.format(i) for i in subList]
                outList.extend(subList)
            else:
                out = self._getValueStr(k, v)
                # out = '> {}'.format(out)
                out = '{}'.format(out)
                outList.append(out)
        finalOut = '\n'.join(outList)
        return finalOut

    def bannerStr(self):
        lines = []
        previous = None
        for k, v in self.items():
            if previous == 'DotMap':
                lines.append('-')
            out = ''
            if isinstance(v, DotMap):
                name = k
                subMap = v
                out = self._getSubMapStr(name, subMap)
                lines.append(out)
                previous = 'DotMap'
            else:
                out = self._getValueStr(k, v)
                lines.append(out)
                previous = 'other'
        lines.append('--')
        s = '\n'.join(lines)
        return s
