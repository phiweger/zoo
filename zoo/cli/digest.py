import click


@click.command()
def encode():
    '''
    \b
    options:
    -- tree
    -- sam
    -- msa
    -- secondary
    overload this function and dispatch to call function from other file
    '''
    pass


@click.command()
def decode():
    '''
    Reconstruct a phylogenetic tree from the database.
    '''
    pass
