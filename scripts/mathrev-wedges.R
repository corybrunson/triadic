# Global diagnostics of Mathematical Reviews network,
# within the aggregate, pure, and applied subnetworks,
# along a fixed-duration sliding window
library(bitriad)
source('code/triadic-spec.R')
source('code/mathrev2igraph.R')
load('calc/mathrev.RData')

#mathrev <- mathrev[sample(nrow(mathrev), 10000), ] # TRIAL RUNS

pure2 <- sprintf('%02d', 3:59)
applied2 <- sprintf('%02d', 60:96)

wedges.df <- function(bigraph) {
    graph <- actor.projection(bigraph)
    w.s.wedges <- choose(degree(graph), 2)
    w.s.closed <- transitivity(graph, type = 'local') * w.s.wedges
    w.s.closed[which(is.na(w.s.closed))] <- 0
    df1 <- data.frame(Wedges = w.s.wedges, Closed = w.s.closed)
    df1$Closed[is.nan(df1$Closed)] <- 0
    df2 <- as.data.frame(do.call(cbind, lapply(
        c(
            injequ.wedges, # Opsahl
            indstr.wedges, # Exclusive
            injstr.wedges,
            injact.wedges,
            indequ.wedges
         ),
        function(f) {
            an.transitivity(bigraph, type = 'both', wedges.fn = f)
        })))
    df <- cbind(df1, df2)
    names(df) <- paste0(rep(c('Classical',
                              'Opsahl',
                              'Exclusive',
                              'InjectiveStructural',
                              'InjectiveActor',
                              'InducedEqual'), each = 2),
                        '.', names(df))
    return(df)
}

mathrev.wedges <- lapply((ran[1] + dur - 1):ran[2], function(yr) {
    wedges.df(as.an(paper.author.graph(mathrev[
        mathrev$year %in% (yr - dur + 1):yr, ])))
})
print('Aggregate done')

mathrev.pure <- mathrev[substr(mathrev$pclass, 1, 2) %in% pure2, ]
mathrev.pure.wedges <- lapply((ran[1] + dur - 1):ran[2], function(yr) {
    wedges.df(as.an(paper.author.graph(mathrev.pure[
        mathrev.pure$year %in% (yr - dur + 1):yr, ])))
})
print('Pure done')

mathrev.applied <- mathrev[substr(mathrev$pclass, 1, 2) %in% applied2, ]
mathrev.applied.wedges <- lapply((ran[1] + dur - 1):ran[2], function(yr) {
    wedges.df(as.an(paper.author.graph(mathrev.applied[
        mathrev.applied$year %in% (yr - dur + 1):yr, ])))
})
print('Applied done')

save(mathrev.wedges, mathrev.pure.wedges, mathrev.applied.wedges,
     file = 'calc/mathrev-wedges.RData')

rm(list = ls())
