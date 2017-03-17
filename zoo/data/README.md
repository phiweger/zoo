## Example data

This data is meant to demonstrate some of the functionality of zoo. Albeit small, we tried to find heterogeneous and somewhat messy data to approximate real use cases.

Access to the package data:

```

get_schema('example.fa')
# /returns/system/path/to/example.fa

from zoo import get_data
from pyfaidx import Fasta

fa = Fasta(get_data('vanilla.fa'))
```

### Influenza A virus, fasta files

[data source](ftp://ftp.ncbi.nih.gov/genomes/INFLUENZA/)

### Deformed wing virus, assembly graph 

- source
- use case

### Ebola

TODO
