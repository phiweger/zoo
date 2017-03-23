![](https://github.com/viehwegerlib/zoo/blob/master/zoo/img/drops.png)

## zoo

A distributable viral database for rapid prototyping (still pre-alpha release).

### Install

```
pip install git+https://github.com/viehwegerlib/zoo.git@master
pip uninstall zoo
```

### Documentation

zoo provides a command-line tool as well as a Python library. For detailed information about zoo's intention, implementation, use cases and tutorials refer to the [wiki](https://github.com/viehwegerlib/zoo/wiki). Documentation about zoo's functions and API is available at [readthedocs](https://readthedocs.org/).

### Quick example

```
zoo import \
    --dat "somelink7234hd8..." \
    --dbname "mock" \
    --client "localhost:27017"

zoo sample \
    --stratified "host" \
    --train-test 30 --split 0.7 \
    --fields "sequence" \
    --fasta "fields,in,header" \
    -o example.fa

zoo msa -i example.fa -o example.mafft.fa
zoo encode --one-hot-encoding example.mafft.fa -o example.mat
# go fire up sklearn

# Use sourmash library to search minhash signatures in a
# sequence Bloom tree.
zoo minhash -k 16,31 -n 200,400
zoo sbt_index --collection "flavivirus,coronavirus" -o sbt_prefix
zoo sbt_search --fastq metagenome.fastq --sbt sbt_prefix -o besthit.csv

# Use ...
```

### Tests

zoo adheres to pytest's package integration [guidance](http://doc.pytest.org/en/latest/goodpractices.html).

```
# cd into package directory
python setup.py test
```

### License

BSD-3-Clause

Copyright (c) 2016 Adrian Viehweger

[![DOI](https://zenodo.org/badge/84596868.svg)](https://zenodo.org/badge/latestdoi/84596868)



