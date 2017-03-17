import os


# specify package data folder location, stackoverflow, 4519127
_ROOT = os.path.abspath(os.path.dirname(__file__))


def get_data(path):
    '''
    import zoo

    print(zoo.get_data('schema/base.json'))
    # prints (system) path to file
    '''
    return os.path.join(_ROOT, 'data', path)
