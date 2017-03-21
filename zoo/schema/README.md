## Schema

zoo is a sequence centric data structure, which influences the database schema design. A sequence record comprises 4 main components or "fields": 

- sequence
- metadata
- relative
- derivative

Sequence information is at the center of zoo's functionalities. It is defined as a string of an arbitrary alphabet, typically RNA, DNA or Protein. As a consequence, each genome segment of a segmented virus such as Influenza A receives its own document, and is linked to the other segments of a given sample.  
Metadata describe the way a sequence "came to be known". Where was it sampled from, who by, from which host, through which sample preparation and sequencing methods?  
Relative information includes taxonomy, phylogeny and linked information. It addressed the question of how a given sequence string compares to others. Parts of phylogenetic trees or multiple sequence alignments are archived in this category.  
Derived information summarizes or reexpresses the information contained in the sequence, including annotations, minhashes and alternative encodings. Derived information is usually heavily dependent on the original sequence. For example, the annotation open reading frame (ORF) derives from the sequence's start and a stop codon position. By definition, derived sequence information  does not by itself make any sense without the underlying raw sequence information.  
Note that all categories interact, e.g. we could use minhash signatures (derivative) to compare a sequence to other ones in the database, storing the top 5 closest sequence ids (relative).

## Workink with schemas

Since zoo is not restricted to any particular use case, some schemas are presented here to give an idea how data can be structured in zoo for different viruses and sets of sequences.

Note that the base schema is "assumed" by most of the zoo package functions, so if you tinker with it (especially deleting keys) things might break. Which is ok, you can probably fix them.

Access to the schemas:

```
from zoo import get_schema
with open(get_schema('base.json')) as infile:
    schema = json.load(infile)
```

## Technicalities

- representing "empty" in JSON, tl;dr: there is no standard (stackoverflow, 21120999)
- ...