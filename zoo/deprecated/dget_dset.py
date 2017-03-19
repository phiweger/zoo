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