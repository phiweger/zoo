## Example data

This data is meant to demonstrate some of the functionality of zoo. Albeit small, we tried to find heterogeneous and somewhat messy data to approximate real use cases.

Access to the package data:

```
# stackoverflow, 779495
import pkg_resources
DATA_PATH = pkg_resources.resource_filename('zoo', 'data/')
# do stuff
```

### Influenza A virus, fasta files

[data source](ftp://ftp.ncbi.nih.gov/genomes/INFLUENZA/)

### Deformed wing virus, assembly graph 

- source
- use case