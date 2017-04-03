import click


@click.command()
def load():
    '''
    Example:

    \b
    zoo load --ncbi accession.txt result.json
    '''
    print('Load.')


@click.command()
def dump():
    '''
    Example:

    \b
    zoo dump --db x --cell y --fasta dump.fa  # all metadata in header, | delim
    '''
    print('Dump.')
