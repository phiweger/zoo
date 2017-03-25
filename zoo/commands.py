import click


@click.option('--count', default=3, help='File to import.')
@click.option('--file', default='foo.bar', help='woof.')
def load(count, file):
    for i in range(count):
        click.echo('Hello %s!' % file)
