import functools
import json


class SeqDoc:
    '''
    - stackoverflow, 5290715, this implementation mixes some answers
    - for object from json: stackoverflow, 6578986

    Example:

    from pymongo import MongoClient
    from zoo.objects import SeqDoc

    client = MongoClient("localhost:27017")
    db = client["zoo"]
    collection = db.get_collection('influenza_a_virus')

    doc = SeqDoc(**collection.find_one())
    doc.metadata.country
    # 'Neverland'
    '''
    def __init__(self, **hashmap):
        '''
        Recursive initialization from a dict.

        Example:

        d = {'a': {'b': [1, 2, 3], 'c': {
            'd': 4}, 'd': [{'hello': 'A'}, {'hello': 'B'}]}, 'b': {}}
        doc = SeqDoc(**d)
        '''
        for key, value in hashmap.items():
            print(key, value, bool(value))
            if isinstance(value, dict) and value != {}:
                # second condition because otherwise json serialization
                # self.json() breaks w/ empty dicts being initialized
                # to SeqDoc object
                self.__dict__[key] = SeqDoc(**value)
            else:
                self.__dict__[key] = value

    def dget_attr(self, attr):  # d .. document
        '''
        modelled after example: stackoverflow, 31174295, much reduced

        Example:

        doc.dget_attr('a.d')  # same as doc.a.d
        # [{'hello': 'A'}, {'hello': 'B'}]
        [i['hello'] for i in doc.dget_attr('a.d')]
        # ['A', 'B']
        '''
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

    def __repr__(self):
        '''
        json.dumps(d, indent=True)
        # return self._id, [i['name'] for i in self.annotation]

        id, names of annotations, seq start, available derivatives, ...
        minhash (k, n tuples)
        msa ids
        '''
        if isinstance(self, SeqDoc) and not isinstance(self, dict):
            return self.__repr__()
        else:
            return self.__dict__



        try:
            # w/ fields of base schema
            pass
        except AttributeError:
            # pretty print complete record in json (=default?)
            # use self.to_json()
            pass

    def __str__():
        '''
        try print seq, else fall back to __repr__
        '''
        pass

    def to_json():
        '''
        We need not implement from_json(), because this (via a hashmap) is
        already done to initialize a SeqDoc object.

        # stackoverflow, 3768895, 10252010
        '''
        pass

    def to_fasta():
        pass

    def get_attr():
        '''
        sd.get_attr('a.b.c')  # gets attribute
        sd.get_attr('a.b.c', 'name')  # w/ annotations

        d = {'a': {'b': [1, 2, 3], 'c': {
            'd': 4}, 'd': [{'hello': 'A'}, {'hello': 'B'}]}}
        sd = SeqDoc(**d)
        sd.a.d
        # [{'hello': 'A'}, {'hello': 'B'}]
        [i['hello'] for i in sd.a.d]
        # ['A', 'B']

        # stackoverflow, 31174295
        rgetattr(sd, 'a.d')
        rsetattr(sd, 'a.b', [2,3,4])
        [i['hello'] for i in rgetattr(sd, 'a.d')]
        '''
        pass


    def to_fasta():
        '''
        return string
        ">la|le|lu\nACTAAGGT..."

        Maybe use a json mask of some sort, would be elegant, i.e.
        pass a nested json with desired keys and done.

        d = {'a': {'b': [1, 2, 3], 'c': {'d', 4}}}
        sd = SeqDoc(**d)
        # do sth. like getattr(sd, 'a.b') - this does not work
        "mask_json"
        ideal:
        header=['metadata.country', 'somefixedlabel', metadata.geo.long]
        getattr(getattr(sd, 'a'), 'b')
        '''
        pass

    def mask_seq():
        '''
        apply an annotation or user defined interval(s)
        '''
        pass

    def seq_apply(somefunc):
        '''
        takes the sequence and any function (e.g. decode gaps)
        '''
        pass
