# Setup
library(bitriad)
library(ggplot2)
library(grid)
library(reshape2)
library(data.table)
source('code/triadic-base.R')
load('calc/example-census.RData')
load('calc/mathrev-census.RData')

# Tally the pairs (x * y, w + z) of wedge tie strength and number of weak ties,
# from a census
stc.count <- function(ftc, max.xy, max.wz) {
    
    # Upper bounds on x * y and w + z
    if(missing(max.xy)) {
        max.xy <- max(sapply(1:nrow(ftc), function(i) {
            lambda <- indexPartition(i - 1)
            lambda[1] * lambda[2]
        }))
    }
    if(missing(max.wz)) {
        max.wz <- indexPartition(nrow(ftc) - 1)[1] + (ncol(ftc) - 1)
    }
    
    # Tallies of ordered triples, counts arrayed by (x * y, w + z)
    dat <- data.frame(Strength = rep(0:max.xy, each = max.wz + 1),
                      WeakTies = rep(0:max.wz, times = max.xy + 1),
                      Count = 0)
    for(i in 1:nrow(ftc)) {
        if(all(ftc[i, ] == 0)) next()
        lambda <- indexPartition(i - 1)
        for(j in 1:ncol(ftc)) {
            if(ftc[i, j] == 0) next()
            w <- j - 1
            for(k in 1:3) {
                strength <- lambda[k] * lambda[(k %% 3) + 1]
                weakties <- w + lambda[((k + 1) %% 3) + 1]
                row <- strength * (max.wz + 1) + weakties + 1
                dat[row, 3] <- dat[row, 3] + ftc[i, j]
            }
        }
    }
    dat
    
}

# STC data for DG1 and GWF
wh.ftc <- which(names(example.census) %in% c('DG1', 'GWF'))
example.stc <- do.call(rbind, lapply(wh.ftc, function(i) {
    cbind(Network = names(example.census)[i], stc.count(example.census[[i]]))
}))

# Maximum strength to include in MR analysis
max.strength <- 20
# STC data for 3 MR intervals
mathrev.stc <- do.call(rbind, lapply(1:3, function(i) {
    cbind(Network = paste0('MR (', yrs[i] - 2, '-', substr(yrs[i], 4, 4), ')'),
          stc.count(mathrev.census[[i]]))
}))
# Restrict to maximum strength
mathrev.stc <- mathrev.stc[mathrev.stc$Strength <= max.strength, ]

# Combine into a single data.table! (woo!)
stc.dt <- as.data.table(rbind(example.stc, mathrev.stc))
# Introduce variable for existence of a weak tie
stc.dt$WeakTie <- stc.dt$WeakTies > 0
# Aggregate the counts over the existence variable
setkey(stc.dt, by = Network,Strength,WeakTie)
stc.dt <- stc.dt[, sum(Count), by = 'Network,Strength,WeakTie']
# Collapse counts
stc.dt <- as.data.table(dcast(stc.dt, Network + Strength ~ WeakTie))
setnames(stc.dt, 'FALSE', 'NoWeakTie')
setnames(stc.dt, 'TRUE', 'WeakTie')
# Add columns for total count, proportion, and standard deviation
stc.dt$Total <- stc.dt[, NoWeakTie + WeakTie]
stc.dt$Proportion <- stc.dt[, WeakTie / Total]
stc.dt$StdDev <- stc.dt[, Proportion * (1 - Proportion) / Total]
# Remove rows with zero Total
stc.dt <- subset(stc.dt, Total > 0)

# Insert empty facet at position 4
stc.dt$Network.alt <- factor(stc.dt$Network,
                             levels = c(levels(stc.dt$Network)[1:2],
                                        '',
                                        levels(stc.dt$Network)[3:5]))

# Plot!
stc.plot <- ggplot(data = stc.dt,
                    aes(x = Strength, y = Proportion)) +
    geom_line() +
    facet_wrap(~ Network.alt, ncol = 3, drop = FALSE, scales = 'free') +
    scale_x_continuous(trans = 'sqrt') +
    theme_bw()

# Remove shaded backdrop from this facet
stc.grob <- ggplotGrob(stc.plot)
## which elements to remove
rm.str <- '(panel|strip_t).*3'
## remove empty panels
stc.grob$grobs[grep(rm.str, names(stc.grob$grobs))] <- NULL
## remove them from the layout
stc.grob$layout = stc.grob$layout[-grep(rm.str, (stc.grob$layout)$name), ]

# Plots
img(width = 2 * wid, height = wid * 4 / 4,
    file = paste0(figloc, 'fig-stc', suf))
grid.draw(stc.grob)
dev.off()

rm(list = ls())
