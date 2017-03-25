import click


@click.option('--file', default='foo.bar', help='File to import.')
def load(file):
    print(file)
