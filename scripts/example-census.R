# Triad censuses of example networks
library(bitriad)
load('calc/example.RData')

example.census <- lapply(example, triad_census_an)
save(example.census, file = 'calc/example-census.RData')

rm(list = ls())
