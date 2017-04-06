library(ggplot2)

df <- read.table('makona.csv', header=T, stringsAsFactors=F, sep=',')
names(df) <- c('date', 'country', 'location', 'latitude', 'longitude')
ggplot(df, aes(x=as.Date(date))) + geom_freqpoly()

