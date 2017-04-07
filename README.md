## zoo - Truly viral messaging.

A portable datastructure for rapid prototyping in (viral) bioinformatics (under development).

![](https://github.com/viehwegerlib/zoo/blob/master/zoo/img/cell.png)

### Install

```
# to use
pip install git+https://github.com/viehwegerlib/zoo.git@master
pip uninstall zoo

# to develop
git clone https://github.com/viehwegerlib/zoo
pip install -e zoo
```

### Documentation

zoo provides a command-line tool as well as a Python library. For detailed information about zoo's intention, implementation, use cases and tutorials refer to the [wiki](https://github.com/viehwegerlib/zoo/wiki/Whitepaper). Documentation about zoo's functions and API is available at [readthedocs](https://readthedocs.org/).

### zoo CLI

```shell
# Stream GenBank records (in JSON format) to data cell, and validate schema.
zoo load --source ncbi --fmt json \
--ids accessions.txt --stdout - | \
zoo init --db mockA --cell foo --validate -
# ... Initializing data cell.
# ... 42 entries inserted into cell "original".
# ... Primary key assigned to field "_id".
# ... inspect cell and commit

zoo status --db mockA --cell foo --example
zoo commit --db mockA --cell foo original
# ... Dumping data cell.
# ... | 42 Elapsed Time: 0:00:00
# ... Done.

# share
mkdir send
cp original.json send/
dat share send/
# ... Syncing Dat Archive: .../send
# ... Link: dat://73401e1b931164763eccsomelonglinkcefc718ebf49f6b4fe4dbad7

# In a faraway place, our collaborator (B) clones a copy of our cell and adds
# it to her "zoo" of other data cells.
mkdir receive
dat clone <link> receive/
zoo add --db mockB --cell foo --primkey genbank.accession receive/original.json
# ... Loading data cell.
# ... Index created on field "genbank.accession".
# ... 39 documents inserted in cell "foo".
# ... 3 duplicates skipped.

# Meanwhile, original.json was modified. B want his zoo to reflect the changes:
dat pull receive/

# diff it
zoo diff --db mockB --cell foo --out diff.json receive/modified.json
zoo pull --db mockB --cell foo receive/modified.json
# ... Updating cell's md5 hashes.
# ... / 0 Elapsed Time: 0:00:00
# ... 
# ... 38 entries unchanged.
# ... 4 entries replaced.

# Now put data cells into your favourite analysis workflow, then use zoo's API 
# to import/ export the results, like multiple sequence or reference-based
# alignments, phylogenetic trees, secondary structure ... happy exploratory
# data analysis. Also, set global vars to reduce typing.
ZOODB=mockB
ZOOCELL=foo
zoo digest --encode tree.nexus
zoo digest --decode msa.mafft.fa

# Not yet implemented: Send metadata about cell to a registry, so others can
# discover it.
zoo push ...

# Create a sequence Bloom tree (SBT) from the minhash signatures of a given 
# cell.
zoo sbt_index --db ref --cell virus --ksize 16 --nsketch 1000 virusref
# ... Initialize SBT.
# ... Compute minhash signatures for selected documents.
# ... k-mer size: 16, sketch size: 1000
# ... \ 9158 Elapsed Time: 0:01:45
# ... Save SBT.
# ... Done.

# Use sourmash_lib to query other signatures about the cell's SBT.
sourmash sbt_search --ksize 16 virusref query.fa.sig
# ... running sourmash subcommand: sbt_search
# ... loaded query: survey.fa... (k=16, DNA)
# ... 0.11 0ef85591-d464-4953-915f-f673907b7e8e (here Zika reference genome)

# Done, lets get some coffee.
zoo drop --db mockB --cell foo --force
zoo destroy --db mockB --force
```

### zoo API

zoo includes a Python library with an intuitive API and many functions to move data in and out of a data cell as well as for formatting the data for downstream analysis, such as:

- export to fasta
- one-hot-encode the selected sequences for machine learning

### Sharing data

As illustrated above, [data cells](https://github.com/viehwegerlib/zoo/wiki/Whitepaper) are shared with the [dat protocol](https://github.com/datproject/dat). It couldn't be easier. Let's say you had a file `zika.json` with some experimental Zika data.

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
# cd into package directory and virtualenv (Python 3)
python setup.py test
```

### License

BSD-3-Clause

Copyright (c) 2016 Adrian Viehweger

[![DOI](https://zenodo.org/badge/84596868.svg)](https://zenodo.org/badge/latestdoi/84596868)



