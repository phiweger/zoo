## Ebola use case

- [source](https://github.com/ebov/space-time) Makona data set
- [source](https://github.com/nextstrain/fauna/blob/master/source-data/geo_lat_long.tsv) geolocation data (longitude, latitude etc.)
- some more [info](https://github.com/nextstrain/fauna/blob/master/EBOLA.md) on how this data is used in the (awesome) [nextstrain](http://www.nextstrain.org/) project

Get data, e.g. by:

```
git clone https://github.com/ebov/space-time
```

We chose this use case, because it is a nice example of what zoo tries to do: Enable the rapid distribution of hand-crafted, lovingly curated datasets, without the need to rewrangle and clean the data, remove all inconsistencies etc., but to get to work straight away. From the 

> This project is a collaboration between the many groups that sequenced virus genomes during the 2014-2016 epidemic in West Africa. __Most of the sequence data has been deposited in GenBank but this repository represents a comprehensive data set curated and annotated by the groups involved.__ The data is currently being analysed with the intent to publish a paper on the spatial and temporal dynamics inferred from virus genome sequences. The analyses and the methods and scripts that underly them will be posted to this repository in the spirit of Open Science and we welcome comments. We hope that early access to this data set and the downstream products such as time-calibrated trees and inferred spatial patterns will foster further research in this area. -- [README, Manoka dataset repo](https://github.com/ebov/space-time)