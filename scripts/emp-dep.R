# Setup
library(bitriad)
library(ggplot2)
library(reshape2)
source('code/triadic-base.R')
load('calc/example.RData')
load('calc/mathrev-wedges.RData')
load('calc/mathrev-closes.RData')
load('calc/mathrev-closed.RData')

# Which statistics to include
stats <- c(project.transitivity,
           opsahl.transitivity,
           indequ.transitivity,
           excl.transitivity)
stat.variety <- c('Classical', 'Opsahl', 'InducedEqual', 'Exclusive')
ex.incl <- c(1, 2, 4)
mr.incl <- c(1, 2, 4)
incl.dyn <- FALSE

# Single data frame of wedge-dependent local triadic closure
example.incl <- which(names(example) %in% c('DG1', 'GWF'))
example.wedges <- do.call(rbind, unlist(lapply(example.incl, function(i) {
    lapply(ex.incl, function(j) {
        wedges <- stats[[j]](example[[i]], type = '')
        data.frame(Network = names(example)[i],
                   Actor = V(example[[i]])$name[!V(example[[i]])$type],
                   Variety = stat.variety[j],
                   Wedges = 2 * wedges[, 1],
                   Alcoves = 2 * wedges[, 2],
                   Closure = wedges[, 2] / wedges[, 1])
    })
}), recursive = FALSE))

# Scatterplot of closure rate versus number of wedges
dep.plots <- ggplot(data = example.wedges,
                    aes(x = Wedges, y = Closure,
                        label = Actor)) +
    geom_point() +# geom_text(size = 3) +
    facet_wrap( ~ Network + Variety, scales = 'free', ncol = 3)

# Plots
img(width = 2 * wid, height = wid * 4 / 4,
    file = paste0(figloc, 'fig-ex-dep', suf))
print(dep.plots + theme_bw())
dev.off()

# Maximum wedge count for MR
l.max <- triangular(7)

# Wedge dependency table
wedge.dependency <-
    function(wedges) {
        agg <- aggregate(wedges[, 2], by = list(wedges[, 1]),
                         FUN = function(x) c(length(x), mean(x)))
        dat <- as.data.frame(cbind(agg$Group.1, agg$x))
        names(dat) <- c('l', 'n', 'C')
        dat$C <- dat$C / dat$l
        dat[which(dat$l > 0), ]
    }

# Make data frames of wedge-dependent local triadic closure
wh.incl <- sapply(stat.variety[mr.incl], function(name) {
    grep(name, names(mathrev.wedges[[1]]))
})
stopifnot(all(dim(wh.incl) == c(2, length(mr.incl))))
dats <- c(
    lapply(1:ncol(wh.incl), function(i) {
        # If 'Classical', recalibrate x-axis
        l.ran <- if(stat.variety[mr.incl][i] == 'Classical') {
            triangular(1:undo.triangular(l.max))
        } else 1:l.max
        # Data frame of wedge counts and average closure rates for each
        data.frame(
            Variety = stat.variety[mr.incl][i],
            Interval = paste0(rep(yrs - 2, each = length(l.ran)), '-',
                         substr(rep(yrs, each = length(l.ran)), 4, 4)),
            Wedges = l.ran,
            Closure = as.vector(sapply(inds, function(j) {
                dat <- wedge.dependency(mathrev.wedges[[j]][, wh.incl[, i]])
                app <- if(stat.variety[mr.incl][i] == 'Classical') {
                    setdiff(triangular(1:undo.triangular(l.max)), dat$l)
                } else {
                    setdiff(1:l.max, dat$l)
                }
                if(length(app) == 0) return(dat$C[which(dat$l <= l.max)])
                sdat <- rbind(dat, data.frame(l = app, n = 0, C = NaN))
                odat <- sdat[order(sdat$l), ]
                return(odat$C[which(odat$l <= l.max)])
            }))
        )
    }),
    if(!incl.dyn) {
        list()
    } else {
        list(data.frame(
            Variety = 'Dynamic',
            Interval = paste0(rep(yrs - 2, each = l.max), '-',
                              substr(rep(yrs, each = l.max), 4, 4)),
            Wedges = 1:l.max,
            Closure = as.vector(sapply(inds, function(j) {
                dat <- wedge.dependency(mathrev.closed[[j]])
                app <- setdiff(1:l.max, dat$l)
                if(length(app) == 0) return(dat$C[which(dat$l <= l.max)])
                sdat <- rbind(dat, data.frame(l = app, n = 0, C = NaN))
                odat <- sdat[order(sdat$l), ]
                return(odat$C[which(odat$l <= l.max)])
            }))
        ))
    }
)

# Conjoin data frames
dat <- do.call(rbind, dats)

# Double wedge counts (since implemented stats take quotients by symmetry)
dat$Wedges <- 2 * dat$Wedges

# Plot closure against wedges for each year, by statistic
dep.plot <- ggplot(data = dat,
                   aes(x = Wedges, y = Closure)) +
    geom_line() +
    facet_grid(Variety ~ Interval, scales = 'free_y')

# Plots
img(width = 2 * wid, height = 1.25 * wid * 4 / 4,
    file = paste0(figloc, 'fig-agg-dep', suf))
print(dep.plot + theme_bw())
dev.off()

rm(list = ls())
