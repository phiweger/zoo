## Visualisation

The examples in this folder should get you started on the most common visualisation tasks.

zoo is not designed to be an end-to-end data-to-insight application. However, its API aims to make data analysis quick and easy, including visualisation (vis).

For us, the [ggplot2](http://docs.ggplot2.org/current/) library of the [R language](https://www.r-project.org/) offers the best data visualisation language, ("gg" stands for "grammar of graphics"). ggplot2 is designed around the [tidy data](https://www.jstatsoft.org/article/view/v059i10) concept. For example, transforming "wide" to "long" data to a flat file in `.csv` format makes life easy. zoo's export API keeps this in mind. 

It is arguable, whether an interactive vis provides more insight than the static images ggplot2 provides. Certainly, with devices like [small multiples](https://en.wikipedia.org/wiki/Small_multiple) dynamic data can be displayed equally well if not more effective.

There are awesome extensions that expand ggplot and its grammar to

- phylogeny: [ggtree](https://github.com/GuangchuangYu/ggtree)
- maps: [ggmap](https://github.com/dkahle/ggmap)
- genomics: [ggbio](http://www.tengfei.name/ggbio/)
- networks: [ggnet2](https://briatte.github.io/ggnet/)

Much of the visual aspect of the vis has been influences by the [bayesplot](https://cran.r-project.org/web/packages/bayesplot/index.html) package and by insights from [E. Tufte](https://en.wikipedia.org/wiki/Edward_Tufte).

