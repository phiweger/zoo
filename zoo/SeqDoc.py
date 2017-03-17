class SeqDoc:
    '''
    - stackoverflow, 5290715, this implementation mixes some answers
    - for object from json: stackoverflow, 6578986

    Example:

    client = MongoClient("localhost:27017")
    db = client["zoo"]
    collection = db.get_collection('influenza_a_virus')

    sd = SeqDoc(**collection.find_one())
    sd.metadata.country
    # 'Neverland'
    '''
    def __init__(self, **hashmap):
        '''
        Recursive initialization from a dict.
        '''
        for key, value in hashmap.items():
            if isinstance(value, dict):
                self.__dict__[key] = SeqDoc(**value)
            else:
                self.__dict__[key] = value


