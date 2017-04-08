import click


@click.command()
def sample():
    '''
    automatic sampling
    '''
    print('sample')


@click.command()
def query():
    '''
    pass query and return cursor? (is this even possible?)
    '''
    print('query')
