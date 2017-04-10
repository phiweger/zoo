library(ape)
library(dplyr)
library(ggmap)
library(ggplot2)
library(ggtree)
source('../zoo/vis/theme_default.R')


# Map


df <- read.table('makona.csv', header=T, stringsAsFactors=F, sep=',')
names(df) <- c('date', 'country', 'location', 'latitude', 'longitude')
df$date <- as.Date(df$date)


df <- df %>% mutate(
    month = format(date, "%m"),
    year = format(date, "%Y")
    )
df <- within(df, m <- paste(year, month, sep='-'))


q <- 
    qmplot(
        longitude, latitude, data=df, maptype='toner-lite', 
        geom='blank') +
        geom_point(color='red', size=.3, alpha=.3) +
        facet_wrap(~ m) +
        theme_default() +
        coord_fixed()

fp <- '~/tmp/map.png'
ggsave(filename=fp, plot=q, height=17, width=17, units='cm', dpi=600)


# Tree

tree <- read.nexus('Makona_1610_genomes_2016-06-23.ml.tree')
# from zoo/data/ebola
p <- ggtree(tree)
fp <- '~/tmp/tree.png'
ggsave(filename=fp, plot=p, height=17, width=17, units='cm', dpi=600)