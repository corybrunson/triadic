# Setup
library(bitriad)
library(ggplot2)
library(grid)
source('code/triadic-base.R')
load('calc/example-census.RData')
load('calc/mathrev-census.RData')

# Tally the pairs (x * y, w + z) of wedge tie strength and number of weak ties,
# from a census
stc.count <- function(ftc, max.xy, max.wz) {
  
  # Upper bounds on x * y and w + z
  if(missing(max.xy)) {
    max.xy <- max(sapply(1:nrow(ftc), function(i) {
      lambda <- index_partition(i - 1)
      lambda[1] * lambda[2]
    }))
  }
  if(missing(max.wz)) {
    max.wz <- index_partition(nrow(ftc) - 1)[1] + (ncol(ftc) - 1)
  }
  
  # Tallies of ordered triples, counts arrayed by (x * y, w + z)
  dat <- data.frame(Strength = rep(0:max.xy, each = max.wz + 1),
                    WeakTies = rep(0:max.wz, times = max.xy + 1),
                    Count = 0)
  for(i in 1:nrow(ftc)) {
    if(all(ftc[i, ] == 0)) next()
    lambda <- index_partition(i - 1)
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

# Combine into a single data frame
stc.df <- dplyr::bind_rows(example.stc, mathrev.stc)
# Introduce variable for existence of a weak tie
stc.df$WeakTie <- stc.df$WeakTies > 0
# Aggregate the counts over the existence variable
stc.df <- dplyr::summarise(dplyr::group_by(stc.df, Network, Strength, WeakTie),
                           Count = sum(Count))
# Cast by existence of weak tie
stc.df <- reshape2::dcast(stc.df,
                          Network + Strength ~ WeakTie,
                          value.var = "Count")
names(stc.df)[3:4] <- c("NoWeakTie", "WeakTie")
# Add columns for total count, proportion, and standard deviation
stc.df <- dplyr::mutate(stc.df, Total = NoWeakTie + WeakTie)
stc.df <- dplyr::mutate(stc.df, Proportion = WeakTie / Total)
stc.df <- dplyr::mutate(stc.df, StdDev = Proportion * (1 - Proportion) / Total)
# Remove rows with zero Total
stc.df <- dplyr::filter(stc.df, Total > 0)

# Insert empty facet at position 3
stc.df$Network.alt <- factor(stc.df$Network,
                             levels = c(unique(stc.df$Network)[1:2],
                                        '',
                                        unique(stc.df$Network)[3:5]))

# Prevent decimal tick marks
brksfn <- function(x) {
  round(waiver(x))
}
# Plot!
stc.plot <- ggplot(data = stc.df,
                   aes(x = Strength, y = Proportion)) +
  geom_line() +
  facet_wrap(~ Network.alt, ncol = 3, drop = FALSE, scales = 'free') +
  scale_x_continuous(trans = 'sqrt', breaks = (1:4) ^ 2) +
  theme_bw()

# Manipulate graphics object table
stc.table <- ggplot_gtable(ggplot_build(stc.plot))
## remove elements (trial and error)
wh_rm <- c(6, 19, 20, 37)
stc.table$grobs <- stc.table$grobs[-wh_rm]
stc.table$layout <- stc.table$layout[-wh_rm, ]

# Plots
img(width = 2 * wid, height = wid * 4 / 4,
    file = paste0(figloc, 'fig-stc', suf))
grid.draw(stc.table)
dev.off()

rm(list = ls())
