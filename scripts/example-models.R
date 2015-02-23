# Transitivity of example networks rel. to null bipartite configuration model
library(bitriad)
library(networksis)
source('code/random.bipartite.R')
load('calc/example.RData')
####
# INCORPORATE THE GLOBAL CALCULATIONS INTO THE LOCAL LOOP;
# SAVE SEPARATELY BUT CONSTRUCT TOGETHER
####

N <- 1000  # total number of simulations
n <- 100   # number of simulations at a time
# Which example networks to compute diagnostics for
# DDGGS2: Example
# DDGGS1: Compare to Opsahl and Liebig-Rao
# BB: Compare to Tutzauer
# GWF: Compare to Faust (1997)
# NMT1: Compare to Liebig-Rao
# NMT2: Include in figure
# FH: Compare to Han
example.incl <- which(names(example) %in% c(
    'DDGGS2'
    , 'DDGGS1'
    #, 'BB'
    #, 'NMT1'
    #, 'NMT2'
    , 'GWF'
    #, 'FH'
))

# Use only first actor in each structural equivalence class
example.reps <- lapply(example, local.reps)
names(example.reps) <- names(example)

# SHOULD HAVE MADE OBSERVATION DATA FRAME SEPARATE FROM LIST OF
# NULL SAMPLE DATA FRAMES

# Local classical, Opsahl, and exclusive clustering
example.model <- list()
for(i in example.incl) {
    print(names(example)[i])
    graph <- actor.projection(example[[i]])
    #vc <- vcount(graph)
    vs1 <- example.reps[[i]]
    vs2 <- which(V(example[[i]])$type == 0)[vs1] # CAUTION
    # Compute local diagnostics on empirical networks
    C.vec <- c(transitivity(graph, type = 'global'),
               transitivity(graph, type = 'local', vids = vs1))
    C.O.vec <- c(opsahl.transitivity(example[[i]], type = 'global'),
                 opsahl.transitivity(example[[i]], type = 'local', vids = vs2))
    C.X.vec <- c(excl.transitivity(example[[i]], type = 'global'),
                 excl.transitivity(example[[i]], type = 'local', vids = vs2))
    # Initialization of data frames, one for each actor
    dfs <- lapply(0:length(vs1), function(v) {
        data.frame(C = C.vec[v + 1],
                   C.O = C.O.vec[v + 1],
                   C.X = C.X.vec[v + 1],
                   W = NA)
    })
    while(nrow(dfs[[1]]) <= N) {
        # Run n simulations at a time
        sim <- igraph.sis(example[[i]], nsim = n)
        # Compute C, C.O, and C.X
        C.mat <- sapply(sim[[1]], function(s) {
            s.proj <- actor.projection(s)
            c(transitivity(s.proj, type = 'global'),
              transitivity(s.proj, vids = vs1, type = 'local'))
        })
        C.O.mat <- sapply(sim[[1]], function(s) {
            wedges <- opsahl.transitivity(s, vids = vs2, type = '')
            c(sum(wedges$T) / sum(wedges$V), wedges$T / wedges$V)
        })
        C.X.mat <- sapply(sim[[1]], function(s) {
            wedges <- excl.transitivity(s, vids = vs2, type = '')
            c(sum(wedges$T) / sum(wedges$V), wedges$T / wedges$V)
        })
        S = 1 / exp(sim$log.prob)
        # Combine with importance wts
        dfs.app <- lapply(0:length(vs1), function(v) {
            data.frame(C = C.mat[v + 1, ],
                       C.O = C.O.mat[v + 1, ],
                       C.X = C.X.mat[v + 1, ],
                       W = S / sum(S))
        })
        dfs <- lapply(1:length(dfs), function(j) rbind(dfs[[j]], dfs.app[[j]]))
        names(dfs) <- c('global', V(graph)$name[vs1])
        print(nrow(dfs[[1]]) - 1)
    }
    example.model <- c(example.model, list(dfs))
    save(example.model, file = 'calc/example-model.RData')
}
names(example.model) <- names(example)[example.incl]
save(example.model, file = 'calc/example-model.RData')

for(pkg in c('networksis', 'ergm', 'network')) {
    detach(paste0('package:', pkg), unload = TRUE)
}
rm(list = ls())
