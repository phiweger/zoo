## Schema

Since zoo is not restricted to any particular use case, some schemas are presented here to give an idea how data can be structured in zoo for different viruses and sets of sequences.

Note that the base schema is "assumed" by most of the zoo package functions, so if you tinker with it (especially deleting keys) things might break. Which is ok, you can probably fix them.


Access to the schemas:

```
# stackoverflow, 779495
import pkg_resources
SCHEMA_PATH = pkg_resources.resource_filename('zoo', 'schema/')
with open(SCHEMA_PATH + 'base.json') as infile:
    schema = json.load(infile)
```

