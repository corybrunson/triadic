# Triad censuses of Mathematical Reviews network,
# within the aggregate, pure, and applied subnetworks,
# over evenly-spaced intervals of equal duration
library(igraph)
source('code/mathrev2igraph.R')
source('code/triadic-spec.R')
load('calc/mathrev.RData')

pure2 <- sprintf('%02d', 3:59)
applied2 <- sprintf('%02d', 60:96)

# Initialize list
mathrev.degree <- list()
mathrev.degree[[1]] <- lapply(years[-(1:2)], function(yr) {
    bigraph <-
        paper.author.graph(mathrev[mathrev$year %in% (yr - dur + 1):yr, ])
    list(actor = tabulate(degree(bigraph)[which(!V(bigraph)$type)]),
         event = tabulate(degree(bigraph)[which(V(bigraph)$type)]))
})

mathrev.pure <- mathrev[substr(mathrev$pclass, 1, 2) %in% pure2, ]
mathrev.degree[[2]] <- lapply(years[-(1:2)], function(yr) {
    bigraph <- paper.author.graph(mathrev.pure[
        mathrev.pure$year %in% (yr - dur + 1):yr,
        ])
    list(actor = tabulate(degree(bigraph)[which(!V(bigraph)$type)]),
         event = tabulate(degree(bigraph)[which(V(bigraph)$type)]))
})

mathrev.applied <- mathrev[substr(mathrev$pclass, 1, 2) %in% applied2, ]
mathrev.degree[[3]] <- lapply(years[-(1:2)], function(yr) {
    bigraph <- paper.author.graph(mathrev.applied[
        mathrev.applied$year %in% (yr - dur + 1):yr,
        ])
    list(actor = tabulate(degree(bigraph)[which(!V(bigraph)$type)]),
         event = tabulate(degree(bigraph)[which(V(bigraph)$type)]))
})

save(mathrev.degree, file = 'calc/mathrev-degree.RData')

rm(list = ls())
