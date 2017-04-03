import click


@click.option('--encode', 'code', flag_value='encode',
              default=True)
@click.option('--decode', 'code', flag_value='decode')
@click.command()
def digest(code):
    '''
    \b
    options:
    -- tree
    -- sam
    -- msa
    -- secondary
    overload this function and dispatch to call function from other file

    \b
    target interface:
    zoo digest --encode tree.nexus
    zoo digest --decode msa.mafft.fa --cell x --id superalignment2017
    '''
    def encode():
        print('We encode.')

    def decode():
        print('We decode.')

    if code == 'encode':
        encode()
    else:
        decode()
