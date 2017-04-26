import click
import json
from zoo import get_schema
from zoo.schema import combine_schema


@click.option('--fp', type=click.Path('r+'), default=None, required=False)
@click.argument('tree')
@click.command()
def schema(fp, tree):
    '''Combine multiple schemas into one validation template.

    Example:

    \b
    zoo schema (a(b,c))  # defaults to zoo's schema directory
    zoo schema --fp path/to/file (a(b,c)) > schema.json
    '''
    if fp is None:
        fp = get_schema()

    print(
        json.dumps(combine_schema(tree, fp), indent=2, sort_keys=True)
        )
