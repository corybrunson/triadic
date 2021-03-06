# Triad censuses of Mathematical Reviews network,
# within the aggregate, pure, and applied subnetworks,
# over evenly-spaced intervals of equal duration
library(bitriad)
source('code/mathrev2igraph.R')
source('code/triadic-spec.R')
load('calc/mathrev.RData')

#mathrev <- mathrev[sample(nrow(mathrev), 1000), ] # TRIAL RUNS

pure2 <- sprintf('%02d', 3:59)
applied2 <- sprintf('%02d', 60:96)

mathrev.applied <- mathrev[substr(mathrev$pclass, 1, 2) %in% applied2, ]
mathrev.applied.census <- lapply(yrs, function(yr) {
  triad_census_an(as_an(paper.author.graph(mathrev.applied[
    mathrev.applied$year %in% (yr - dur + 1):yr, ])))
})
mathrev.applied.unif.census <- lapply(mathrev.applied.census, function(census) {
  project_census(census)$difference
})
save(mathrev.applied.census, mathrev.applied.unif.census,
     file = 'calc/mathrev-applied-census.RData')
print('Applied done')

mathrev.pure <- mathrev[substr(mathrev$pclass, 1, 2) %in% pure2, ]
mathrev.pure.census <- lapply(yrs, function(yr) {
  triad_census_an(as_an(paper.author.graph(mathrev.pure[
    mathrev.pure$year %in% (yr - dur + 1):yr, ])))
})
mathrev.pure.unif.census <- lapply(mathrev.pure.census, function(census) {
  project_census(census)$difference
})
save(mathrev.applied.census, mathrev.pure.unif.census,
     file = 'calc/mathrev-pure-census.RData')
print('Pure done')

mathrev.census <- lapply(yrs, function(yr) {
  triad_census_an(as_an(paper.author.graph(mathrev[
    mathrev$year %in% (yr - dur + 1):yr, ])))
})
mathrev.unif.census <- lapply(mathrev.census, function(census) {
  project_census(census)$difference
})
save(mathrev.census, mathrev.unif.census,
     file = 'calc/mathrev-census.RData')
print('Aggregate done')

rm(list = ls())
