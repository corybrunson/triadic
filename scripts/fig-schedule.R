# Figure: An affiliation network and two of its scheduled subgraphs
library(bitriad)
source('code/triadic-base.R')

# Example graph
g <- graph(c(1,6,1,8,1,9,
             2,6,2,7,
             3,11,3,12,3,8,3,9,3,7,3,6,
             4,6,4,7,4,10,
             5,10), directed = FALSE)
V(g)$name <- c('i', 'j', 'k', 'l', 'm', 1:(vcount(g) - 5))
V(g)$type <- 1:vcount(g) > 5

# Layout
set.seed(100)
lfr <- layout.fruchterman.reingold(g)

# Frame for additional nodes
lfr.xran <- range(lfr[, 1])
lfr.yran <- range(lfr[, 2])
frame.mat <- matrix(c(lfr.xran[1] - .25 * diff(lfr.xran),
                      lfr.xran[1] - .25 * diff(lfr.xran),
                      lfr.xran[2] + .25 * diff(lfr.xran),
                      lfr.xran[2] + .25 * diff(lfr.xran),
                      lfr.yran[1] - .25 * diff(lfr.yran),
                      lfr.yran[2] + .25 * diff(lfr.yran),
                      lfr.yran[1] - .25 * diff(lfr.yran),
                      lfr.yran[2] + .25 * diff(lfr.yran)), nc = 2)

# Plot
img(width = wid, height = 2 * wid,
    file = paste0(figloc, 'fig-sched', suf))
plot(add.vertices(g, 4),
     layout = rbind(lfr, frame.mat),
     vertex.label = c(V(g)$name, rep('', 4)),
     vertex.shape = c(g.shape[V(g)$type + 1], rep('circle', 4)),
     vertex.color = c(g.color[V(g)$type + 1], rep('white', 4)),
     vertex.size = c(g.size[V(g)$type + 1], rep(1, 4)) * .7,
     edge.width = edge.wid,
     edge.color = 'black',
     vertex.label.family = 'sans',
     vertex.label.color = 'white')
dev.off()

# Plot schedule of i, j, k, l
rm.nodes <- c(5, 10, 11, 12)
img(width = wid, height = 2 * wid,
    file = paste0(figloc, 'fig-sched1', suf))
plot(add.vertices(delete.vertices(g, rm.nodes), 4),
     layout = rbind(lfr[-rm.nodes, ], frame.mat),
     vertex.label = c(V(g)$name[-rm.nodes], rep('', 4)),
     vertex.shape = c(g.shape[V(g)$type + 1][-rm.nodes], rep('circle', 4)),
     vertex.color = c(g.color[V(g)$type + 1][-rm.nodes], rep('white', 4)),
     vertex.size = c(g.size[V(g)$type + 1][-rm.nodes], rep(1, 4)) * .7,
     edge.width = edge.wid,
     edge.color = 'black',
     vertex.label.family = 'sans',
     vertex.label.color = 'white')
dev.off()

# Plot schedule of k, l, m
rm.nodes <- c(1, 2, 8, 9, 11, 12)
img(width = wid, height = 2 * wid,
    file = paste0(figloc, 'fig-sched2', suf))
plot(add.vertices(delete.vertices(g, rm.nodes), 4),
     layout = rbind(lfr[-rm.nodes, ], frame.mat),
     vertex.label = c(V(g)$name[-rm.nodes], rep('', 4)),
     vertex.shape = c(g.shape[V(g)$type + 1][-rm.nodes], rep('circle', 4)),
     vertex.color = c(g.color[V(g)$type + 1][-rm.nodes], rep('white', 4)),
     vertex.size = c(g.size[V(g)$type + 1][-rm.nodes], rep(1, 4)) * .7,
     edge.width = edge.wid,
     edge.color = 'black',
     vertex.label.family = 'sans',
     vertex.label.color = 'white')
dev.off()

rm(list = ls())
