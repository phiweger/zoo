## Example data

This data is meant to demonstrate some of the functionality of zoo. Albeit small, we tried to find heterogeneous and somewhat messy data to approximate real use cases.

Access to the package data:

```
from zoo import get_data
from pyfaidx import Fasta

fa = Fasta(get_data('vanilla.fa'))
fa = Fasta(get_data('zika/survey.fa'))
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

### Zika

We will pretend we sampled a couple of samples to detect Zika virus.

- data source: BioProject ID [PRJNA344504](http://www.ebi.ac.uk/ena/data/view/PRJNA344504&portal=sequence_update), downloaded manually as `survey.fa`, 2017-03-24
- `survey.txt`: corresponding GenBank accession codes
- `survey.json`: corresponding data cell (= database dump)
- `survey.sig`: minhash signature of all sequences in data cell
- `survey.py`: see how the files were generated

