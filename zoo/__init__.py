import os


# specify package data folder location, stackoverflow, 4519127
_ROOT = os.path.abspath(os.path.dirname(__file__))


def get_data(path=None):
    '''
    import zoo

    print(zoo.get_data('schema/base.json'))
    # prints (system) path to file
    '''
    if path is None:
        return os.path.join(_ROOT, 'data')
    else:
        return os.path.join(_ROOT, 'data', path)


def get_schema(path=None):
    '''
    import zoo

    print(zoo.get_data('schema/base.json'))
    # prints (system) path to file
    '''
    if path is None:
        return os.path.join(_ROOT, 'schema')
    else:
        return os.path.join(_ROOT, 'schema', path)


def get_script(path):
    '''
    import zoo

    print(zoo.get_data('schema/base.json'))
    # prints (system) path to file
    '''
    return os.path.join(_ROOT, 'scripts', path)
