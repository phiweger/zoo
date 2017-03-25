## zoo

A distributable viral database for rapid prototyping (still pre-alpha release).

![](https://github.com/viehwegerlib/zoo/blob/master/zoo/img/cell.png)

### Install

```
pip install git+https://github.com/viehwegerlib/zoo.git@master
pip uninstall zoo
```

### Documentation

zoo provides a command-line tool as well as a Python library. For detailed information about zoo's intention, implementation, use cases and tutorials refer to the [wiki](https://github.com/viehwegerlib/zoo/wiki/Whitepaper). Documentation about zoo's functions and API is available at [readthedocs](https://readthedocs.org/).

### API teaser (under development)

```
zoo load --json ...
zoo load --ncbi ...
zoo load --ebi ...
zoo load --dat ...

zoo dump --fasta file.fasta # dump entire db, default: only UUID exported to fasta header
zoo dump --fields ... --ids ... --json file.json  # select fields and ids
zoo dump --query q.json ...  # or pass a valid mongodb query directly
zoo dump --pipeline pl.json ...  # same goes for pipelines
dat share .  # to share

zoo sample ... --fasta file.fasta  # code for --fasta, --query etc. recycled

# do MSA using Mafft
zoo msa --attach alignment.mafft.fa "derivative.msa"
zoo msa --reconstruct "ID" --out alignment.mafft.fa

# or align against some reference
zoo sam --attach alignment.sam "derivative.msa"  # note the acronym: sam, msa, + 1
zoo sam --reconstruct "ID" --out alignment.sam

# FastML, RAxML refinement
zoo tree --attach tree.phylip
zoo tree --reconstruct "ID" --out tree.phylip

# minhash
# select, sample to fasta
sourmash compute ...
zoo minhash --attach sourmash.sig
zoo minhash --sbt prefix --collections "list,of,collections"
zoo minhash --sbt prefix --query q.json
```

### Sharing data

[Data cell](https://github.com/viehwegerlib/zoo/wiki/Whitepaper) are shared with the [dat protocol](https://github.com/datproject/dat). It couldn't be easier. Let's say you had a file `zika.json` with some experimental Zika data.

```
dat share .../zika_survey/  # contains zika.json
# Syncing Dat Archive: .../zika_survey
# Link: dat://ff92ce30e1ff6ebd75edeb42f04239367243a58b7838f50706bd995e5dbc5d4c
```

We can send this link to a colleague or put it in the zoo registry for others to find.

```
# meanwhile in a faraway place
dat clone ff92ce30e1ff6ebd75edeb42f04239367243a58b7838f50706bd995e5dbc5d4c
ls
# zika.json
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



