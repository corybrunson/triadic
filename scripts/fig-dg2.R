# Figure: DG2 and its one-mode projection
library(igraph)
library(bitriad)
source('code/triadic-base.R')
load('calc/example.RData')

dg <- example$DG2
dg.proj <- actor_projection(dg)

# Layouts for two-mode network and its one-mode projection
set.seed(10)
lay.dg <- layout.fruchterman.reingold(dg, niter = 100)
lay.proj <- lay.dg[1:5, ]
# Plot the graphs
img(width = wid, height = wid, file = paste0(figloc, 'fig-dg2', suf))
plot(dg, layout = lay.dg,
     vertex.color = ifelse(V(dg)$type == 0, g.color[1], g.color[2]),
     vertex.shape = ifelse(V(dg)$type == 0, g.shape[1], g.shape[2]),
     vertex.size = ifelse(V(dg)$type == 0, g.size[1], g.size[2]),
     edge.width = edge.wid, edge.color = 'black',
     vertex.label = c(LETTERS[1:5], 1:5),
     vertex.label.family = 'sans', vertex.label.color = 'white')
dev.off()
img(width = wid, height = wid,
    file = paste0(figloc, 'fig-dg2proj', suf))
plot(dg.proj, layout = lay.proj,
     vertex.color = g.color[1],
     vertex.shape = g.shape[1],
     vertex.size = g.size[1],
     edge.width = 2, edge.color = 'black',
     vertex.label = LETTERS[1:5],
     vertex.label.family = 'sans', vertex.label.color = 'white')
dev.off()

rm(list = ls())
