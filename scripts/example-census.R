# Triad censuses of example networks
library(bitriad)
load('calc/example.RData')

example.census <- lapply(example, an.triad.census)
save(example.census, file = 'calc/example-census.RData')

rm(list = ls())
