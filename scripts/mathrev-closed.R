# Wedge closure in the Mathematical Reviews network,
# within the aggregate, pure, and applied subnetworks,
# aggregated over adjacent 3-year intervals
library(bitriad)
source('code/triadic-spec.R')
source('code/mathrev2igraph.R')
load('calc/mathrev.RData')

#mathrev <- mathrev[sample(nrow(mathrev), 10000), ] # TRIAL RUNS

pure2 <- sprintf('%02d', 3:59)
applied2 <- sprintf('%02d', 60:96)

mathrev.pure <- mathrev[substr(mathrev$pclass, 1, 2) %in% pure2, ]
mathrev.pure.closed <- lapply((ran[1] + dur - 1):ran[2], function(yr) {
    dyn.triadic.closure.bigraph(
        paper.author.graph(
            mathrev.pure[mathrev.pure$year %in% (yr - dur + 1):yr, ]
        ),
        memory = Inf,
        type = 'both'
    )
})
print('Pure done')

mathrev.applied <- mathrev[substr(mathrev$pclass, 1, 2) %in% applied2, ]
mathrev.applied.closed <- lapply((ran[1] + dur - 1):ran[2], function(yr) {
    dyn.triadic.closure.bigraph(
        paper.author.graph(
            mathrev.applied[mathrev.applied$year %in% (yr - dur + 1):yr, ]
        ),
        memory = Inf,
        type = 'both'
    )
})
print('Applied done')

mathrev.closed <- lapply((ran[1] + dur - 1):ran[2], function(yr) {
    dyn.triadic.closure.bigraph(
        paper.author.graph(
            mathrev[mathrev$year %in% (yr - dur + 1):yr, ]
        ),
        memory = Inf,
        type = 'both'
    )
})
print('Aggregate done')

save(mathrev.closed, mathrev.pure.closed, mathrev.applied.closed,
     file = 'calc/mathrev-closed.RData')

rm(list = ls())
