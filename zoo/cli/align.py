import click


@click.command()
def msa_encode():
    '''
    Represent a multiple sequence alignment in the database.
    '''
    pass


@click.command()
def msa_decode():
    '''
    Reconstruct a multiple sequence alignment from the database.
    '''
    pass


@click.command()
def sam_encode():
    '''
    Represent a reference-based alignment in the database.
    '''
    pass


@click.command()
def sam_decode():
    '''
    Reconstruct a reference-based alignment from the database.
    '''
    pass
