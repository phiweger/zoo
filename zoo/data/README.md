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

### Testing minhash and SBT functions

For generation of `mock_flu.fa` see `mock_flu.py.`
 
### Deformed wing virus, assembly graph 

- source
- use case

### Ebola

TODO

### JMTV segmented virus

[From](http://www.pnas.org/content/111/18/6744.full): Qin, X.-C. et al. A tick-borne segmented RNA virus contains genome segments
derived from unsegmented viral ancestors. PNAS 111, 6744–6749 (2014).

See `examples/segmented_virus.py` for data generation and use case.