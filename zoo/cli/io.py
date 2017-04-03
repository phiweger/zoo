import click


@click.command()
def load():
    '''
    fasta, fastq, ...
    '''
    print('Load.')


@click.command()
def dump():
    '''
    fasta, fastq, ...
    '''
    print('Dump.')
