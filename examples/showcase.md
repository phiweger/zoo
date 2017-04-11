Why?

- high level: design, concept
- integration: registry, downstream analysis, visualisation (nextstrain)

```
cp ~/repos/zoo/zoo/data/tests/cell_a.json .
cp ~/repos/zoo/zoo/data/tests/cell_b.json .

echo $ZOODB
echo $ZOOCELL

zoo init --db test --cell jena cell_a.json

# primary key added
zoo status --db test --cell jena --example

zoo add --db test --cell jena cell_b.json

mkdir send
mkdir receive

zoo commit --db test --cell jena --n 5 cell_c
subl cell_c.json
subl cell_c.zoo  # for registry

cp cell_c.json send
dat share send
# cp link

# other window
dat clone link receive
dat pull receive
dat sync receive
```

```
from pymongo import MongoClient
c = MongoClient('localhost:27017')['test']['jena']
```

