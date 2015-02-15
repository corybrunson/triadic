# Figure: The two-mode kite graph and two possible lifts of it
library(bitriad)
source('code/triadic-base.R')

# Graphs: two bipartites that project to the same unipartite
g <- graph(c(1, 2, 2, 3, 3, 1, 3, 4), directed = FALSE)
b1 <- graph(c(1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 1, 6, 7, 7, 8), directed = FALSE)
b2 <- graph(c(1, 2, 1, 3, 1, 4, 4, 5, 5, 6), directed = FALSE)
b3 <- graph(c(1, 2, 2, 3, 3, 4, 4, 5, 4, 6, 5, 7, 6, 7, 7, 1, 7, 8, 8, 9),
            directed = FALSE)
b4 <- graph(c(1, 4, 1, 5, 1, 6, 2, 4, 2, 5, 2, 6, 3, 4, 3, 5, 3, 6, 6, 7, 7, 8),
            directed = FALSE)

# actors (0) and events (1)
V(g)$type <- rep(0, 4)
V(b1)$type <- c(1, 0, 1, 0, 1, 0, 1, 0)
V(b2)$type <- c(1, 0, 0, 0, 1, 0)
V(b3)$type <- c(1, 0, 1, 0, 1, 1, 0, 1, 0)
V(b4)$type <- c(1, 1, 1, 0, 0, 0, 1, 0)

# Label size (wtf)
V(g)$label.cex = .75
V(b1)$label.cex = .75
V(b2)$label.cex = .75
V(b3)$label.cex = .75
V(b4)$label.cex = .75

# Layouts: rotated slightly for aesthetics
l <- matrix(c(
    -cos(0:2 * 2 * pi / 3 + pi / 6), 0,
    sin(0:2 * 2 * pi / 3 + pi / 6), -1 - sqrt(3)
), nc = 2)
l1 <- matrix(c(
    -cos(0:5 * pi / 3 - pi / 6), -.5, 0,
    sin(0:5 * pi / 3 - pi / 6), -1 - sqrt(3) / 2, -1 - sqrt(3)
), nc = 2)
l2 <- matrix(c(
    0, -cos(0:2 * 2 * pi / 3 + pi / 6), -.5, 0,
    0, sin(0:2 * 2 * pi / 3 + pi / 6), -1 - .5 * sqrt(3), -1 - sqrt(3)
), nc = 2)
l3 <- matrix(c(
    -cos(0:3 * pi / 3 - pi / 6),
    -cos(4 * pi / 3 - pi / 6) * (1 - .333),
    -cos(4 * pi / 3 - pi / 6) * (1 + .333),
    -cos(5 * pi / 3 - pi / 6),
    -.5, 0,
    sin(0:3 * pi / 3 - pi / 6),
    sin(4 * pi / 3 - pi / 6) * (1 - .333),
    sin(4 * pi / 3 - pi / 6) * (1 + .333),
    sin(5 * pi / 3 - pi / 6),
    -1 - sqrt(3) / 2, -1 - sqrt(3)
), nc = 2)
l4 <- matrix(c(
    -cos(0:2 * 2 * pi / 3 - pi / 6) * .4,
    -cos(0:2 * 2 * pi / 3 + pi / 6), -.5, 0,
    sin(0:2 * 2 * pi / 3 - pi / 6) * .4,
    sin(0:2 * 2 * pi / 3 + pi / 6), -1 - .5 * sqrt(3), -1 - sqrt(3)
), nc = 2)

# Plots!
img(width = wid, height = 2 * wid,
    file = paste0(figloc, 'fig-kite-proj', suf))
plot(
    add.vertices(g, 4),
    vertex.label = c('i', 'j', 'k', 'l', rep('', 4)),
    vertex.shape = c(g.shape[V(g)$type + 1], rep('circle', 4)),
    vertex.color = c(g.color[V(g)$type + 1], rep('white', 4)),
    vertex.size = c(g.size[V(g)$type + 1], rep(1, 4)) * .5,
    edge.width = edge.wid,
    edge.color = 'black',
    vertex.label.family = 'sans',
    vertex.label.color = 'white',
    asp = .5,
    layout = rbind(l, matrix(c(-4, -4, 4, 4, -3, 1, -3, 1), nc = 2))
)
dev.off()

img(width = wid, height = 2 * wid,
    file = paste0(figloc, 'fig-kite-aff1', suf))
plot(
    add.vertices(b1, 4),
    vertex.label = c(c('', 'i', '', 'j', '', 'k', '', 'l'), rep('', 4)),
    vertex.shape = c(g.shape[V(b1)$type + 1], rep('circle', 4)),
    vertex.color = c(g.color[V(b1)$type + 1], rep('white', 4)),
    vertex.size = c(g.size[V(b1)$type + 1], rep(1, 4)) * .5,
    edge.width = edge.wid,
    edge.color = 'black',
    vertex.label.family = 'sans',
    vertex.label.color = 'white',
    asp = .5,
    layout = rbind(l1, matrix(c(-4, -4, 4, 4, -3, 1, -3, 1), nc = 2))
)
dev.off()

img(width = wid, height = 2 * wid,
    file = paste0(figloc, 'fig-kite-aff2', suf))
plot(
    add.vertices(b2, 4),
    vertex.label = c(c('', 'i', 'j', 'k', '', 'l'), rep('', 4)),
    vertex.shape = c(g.shape[V(b2)$type + 1], rep('circle', 4)),
    vertex.color = c(g.color[V(b2)$type + 1], rep('white', 4)),
    vertex.size = c(g.size[V(b2)$type + 1], rep(1, 4)) * .5,
    edge.width = edge.wid,
    edge.color = 'black',
    vertex.label.family = 'sans',
    vertex.label.color = 'white',
    asp = .5,
    layout = rbind(l2, matrix(c(-4, -4, 4, 4, -3, 1, -3, 1), nc = 2))
)
dev.off()

img(width = wid, height = 2 * wid,
    file = paste0(figloc, 'fig-kite-aff3', suf))
plot(
    add.vertices(b3, 4),
    vertex.label = c(c('', 'i', '', 'j', '', '', 'k', '', 'l'), rep('', 4)),
    vertex.shape = c(g.shape[V(b3)$type + 1], rep('circle', 4)),
    vertex.color = c(g.color[V(b3)$type + 1], rep('white', 4)),
    vertex.size = c(g.size[V(b3)$type + 1], rep(1, 4)) * .5,
    edge.width = edge.wid,
    edge.color = 'black',
    vertex.label.family = 'sans',
    vertex.label.color = 'white',
    asp = .5,
    layout = rbind(l3, matrix(c(-4, -4, 4, 4, -3, 1, -3, 1), nc = 2))
)
dev.off()

img(width = wid, height = 2 * wid,
    file = paste0(figloc, 'fig-kite-aff4', suf))
plot(
    add.vertices(b4, 4),
    vertex.label = c(c('', '', '', 'i', 'j', 'k', '', 'l'), rep('', 4)),
    vertex.shape = c(g.shape[V(b4)$type + 1], rep('circle', 4)),
    vertex.color = c(g.color[V(b4)$type + 1], rep('white', 4)),
    vertex.size = c(g.size[V(b4)$type + 1], rep(1, 4)) * .5,
    edge.width = edge.wid,
    edge.color = 'black',
    vertex.label.family = 'sans',
    vertex.label.color = 'white',
    asp = .5,
    layout = rbind(l4, matrix(c(-4, -4, 4, 4, -3, 1, -3, 1), nc = 2))
)
dev.off()

rm(list = ls())
