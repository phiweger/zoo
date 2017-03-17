import os


# specify package data folder location, stackoverflow, 4519127
_ROOT = os.path.abspath(os.path.dirname(__file__))


def get_data(path):
    '''
    print(get_data('resource1/foo.txt'))
    '''
    return os.path.join(_ROOT, 'data', path)
