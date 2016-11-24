# Wedge closure in the Mathematical Reviews network,
# within the aggregate, pure, and applied subnetworks,
# along a fixed-duration sliding window
library(igraph)
library(bitriad)
source('code/triadic-spec.R')
source('code/mathrev2igraph.R')
load('calc/mathrev.RData')

#mathrev <- mathrev[sample(nrow(mathrev), 10000), ] # TRIAL RUNS

pure2 <- sprintf('%02d', 3:59)
applied2 <- sprintf('%02d', 60:96)

mathrev.closes <- lapply(1:(dur - 1), function(d) {
  lapply((ran[1] + dur - 1):ran[2], function(yr) {
    wedge.closure(mathrev, (yr - dur + 1):(yr - d), (yr - d + 1):yr,
                  type = 'both')
  })
})

mathrev.pure <- mathrev[substr(mathrev$pclass, 1, 2) %in% pure2, ]
mathrev.pure.closes <- lapply(1:(dur - 1), function(d) {
  lapply((ran[1] + dur - 1):ran[2], function(yr) {
    wedge.closure(mathrev.pure, (yr - dur + 1):(yr - d), (yr - d + 1):yr,
                  type = 'both')
  })
})

mathrev.applied <- mathrev[substr(mathrev$pclass, 1, 2) %in% applied2, ]
mathrev.applied.closes <- lapply(1:(dur - 1), function(d) {
  lapply((ran[1] + dur - 1):ran[2], function(yr) {
    wedge.closure(mathrev.applied, (yr - dur + 1):(yr - d), (yr - d + 1):yr,
                  type = 'both')
  })
})

save(mathrev.closes, mathrev.pure.closes, mathrev.applied.closes,
     file = 'calc/mathrev-closes.RData')

rm(list = ls())
