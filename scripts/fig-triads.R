# Figure: Example triads, including canonical open and closed wedges
library(bitriad)
source('code/triadic-base.R')

triad.scale <- 1.2

# Plot example triads
example.triads <- list(
  list(c(0,0,0), 0),
  list(c(0,0,0), 1),
  list(c(0,0,0), 2),
  list(c(0,0,0), 3),
  list(c(1,0,0), 0),
  list(c(1,0,0), 1),
  list(c(1,0,0), 2),
  list(c(1,1,0), 0),
  list(c(1,1,0), 1),
  list(c(1,1,1), 0),
  list(c(1,1,1), 1),
  list(c(2,0,0), 0),
  list(c(2,1,0), 1),
  list(c(2,1,1), 0)
)

for(ex in example.triads) {
  par(mar = rep(0, 4))
  img(height = 4 / triad.scale, width = 3.33 / triad.scale,
      file = paste0(figloc, 'fig-triad-',
                    paste(ex[[1]], collapse = ''), '-', ex[[2]], suf))
  plot_triad(lambda = ex[[1]], w = ex[[2]],
             cex = triad.scale, scale = .25,
             event_names = rep('', sum(unlist(ex))))
  dev.off()
}

for(ex in example.triads) {
  par(mar = rep(0, 4))
  img(height = 4 / triad.scale, width = 3.33 / triad.scale,
      file = paste0(figloc, 'fig-triad-anon-',
                    paste(ex[[1]], collapse = ''), '-', ex[[2]], suf))
  plot_triad(lambda = ex[[1]], w = ex[[2]],
             cex = triad.scale, scale = .25,
             actor_names = rep('', 3),
             event_names = rep('', sum(unlist(ex))))
  dev.off()
}

par(mar = rep(0, 4))
img(height = 4 / triad.scale, width = 3.33 / triad.scale,
    file = paste0(figloc, 'fig-triad-induced1', suf))
plot_triad(lambda = c(1, 0, 0), w = 1,
           cex = triad.scale, scale = .25,
           actors = c('i', 'j', 'k'), events = c('d', 'e'))
dev.off()

par(mar = rep(0, 4))
img(height = 4 / triad.scale, width = 3.33 / triad.scale,
    file = paste0(figloc, 'fig-triad-induced2', suf))
plot_triad(lambda = c(2, 1, 0), w = 1,
           cex = triad.scale, scale = .25,
           actors = c('i', 'j', 'k'), events = c('d', 'f', 'g', 'e'))
dev.off()

# traditional network triads
tr <- make_triad(c(0, 0, 0), 1)
layout <- layout_triad(lambda = c(0, 0, 0), w = 1)
tr2 <- add.edges(delete.vertices(tr, 4), c(1, 2, 2, 3))
tr3 <- add.edges(tr2, c(1, 3))
layout <- layout[1:3, ]
xlim <- c(-1.4, 1.4)
ylim <- c(-1.4, 1.4)

# plot as in plot_triad
par(mar = rep(0, 4))
img(height = 4 / triad.scale, width = 3.33 / triad.scale,
    file = paste0(figloc, 'fig-triad-uni2', suf))
plot(tr2, layout = layout,
     xlim = xlim, ylim = ylim,
     vertex.label = rep('', vcount(tr2)),
     vertex.shape = rep('circle', 3),
     vertex.size = rep(36, 3) * triad.scale,
     vertex.color = rep('SkyBlue2', 3),
     vertex.label.family = 'sans', vertex.label.font = 2,
     vertex.label.color = 'white',
     edge.width = 2, edge.color = 'black', rescale = FALSE, asp = 0)
dev.off()
par(mar = rep(0, 4))
img(height = 4 / triad.scale, width = 3.33 / triad.scale,
    file = paste0(figloc, 'fig-triad-uni3', suf))
plot(tr3, layout = layout,
     xlim = xlim, ylim = ylim,
     vertex.label = rep('', vcount(tr2)),
     vertex.shape = rep('circle', 3),
     vertex.size = rep(36, 3) * triad.scale,
     vertex.color = rep('SkyBlue2', 3),
     vertex.label.family = 'sans', vertex.label.font = 2,
     vertex.label.color = 'white',
     edge.width = 2, edge.color = 'black', rescale = FALSE, asp = 0)
dev.off()

rm(list = ls())
