# Setup
library(bitriad)
library(ggplot2)
library(grid)
library(reshape2)
source('code/triadic-base.R')
source('code/random.bipartite.R')
load('calc/mathrev-wedges.RData')
load('calc/mathrev-closed.RData')
load('calc/mathrev-closes.RData')
load('calc/mathrev-degree.RData')

# Reset years
years <- 1987:2008

# Conjoin wedge proportions for static and dynamic global statistics
global.dat <- do.call(rbind, list(
    data.frame(Type = 'ClusteringCoefficient',
               Variety = 'Classical',
               Year = years,
               Statistic = sapply(mathrev.wedges, wedges2cc, 1, 2))
    , data.frame(Type = 'ClusteringCoefficient',
                 Variety = 'Opsahl',
                 Year = years,
                 Statistic = sapply(mathrev.wedges, wedges2cc, 3, 4))
    , data.frame(Type = 'ClusteringCoefficient',
                 Variety = 'Exclusive',
                 Year = years,
                 Statistic = sapply(mathrev.wedges, wedges2cc, 5, 6))
    #, data.frame(Type = 'ClusteringCoefficient',
    #             Variety = 'InjectiveStructural',
    #             Year = years,
    #             Statistic = sapply(mathrev.wedges, wedges2cc, 7, 8))
    , data.frame(Type = 'ClusteringCoefficient',
                 Variety = 'InjectiveActor',
                 Year = years,
                 Statistic = sapply(mathrev.wedges, wedges2cc, 9, 10))
    , data.frame(Type = 'TriadicClosure',
                 Variety = 'Classical',
                 Year = years,
                 Statistic = sapply(mathrev.wedges, wedges2tr, 1, 2))
    , data.frame(Type = 'TriadicClosure',
                 Variety = 'Exclusive',
                 Year = years,
                 Statistic = sapply(mathrev.wedges, wedges2tr, 5, 6))
    , data.frame(Type = 'TriadicClosure',
                 Variety = 'InjectiveActor',
                 Year = years,
                 Statistic = sapply(mathrev.wedges, wedges2tr, 9, 10))
    , data.frame(Type = 'TriadicClosure',
                 Variety = 'Dynamic',
                 Year = years,
                 Statistic = sapply(mathrev.closed, wedges2tr, 1, 2))
    #, data.frame(Type = 'TriadicClosure',
    #             Variety = 'Dynamic1',
    #             Year = years,
    #             Statistic = sapply(mathrev.closes[[1]], wedges2tr, 1, 2))
    #, data.frame(Type = 'TriadicClosure',
    #             Variety = 'Dynamic2',
    #             Year = years,
    #             Statistic = sapply(mathrev.closes[[2]], wedges2tr, 1, 2))
    ))

# Restrict to adjacent 3-year intervals
global.dat <- global.dat[global.dat$Year %in% seq(1987, 2008, 3), ]

# Percent variation in DTC explained by Exclusive T
ts.dtc <- global.dat$Statistic[global.dat$Variety == 'Dynamic']
ts.tex <- global.dat$Statistic[global.dat$Variety == 'Exclusive' &
                                   global.dat$Type == 'TriadicClosure']
print('Percent of mean value of dynamic TC explained by exclusive TC')
print(round(1 - (sum(ts.dtc) - sum(ts.tex)) / sum(ts.dtc), 2))
print('Percent of variation in dynamic TC explained by exclusive TC')
print(round(1 - var(ts.dtc - ts.tex) / var(ts.dtc), 2))

# Time series plots for global statistics on the aggregate MR network
global.plot <- ggplot(data = global.dat,
                  aes(x = Year, y = Statistic, shape = Variety)) +
    geom_point() + geom_line() +
    facet_wrap(~ Type, scales = 'free') +
    scale_y_log10() + scale_x_continuous('End year') +
    scale_shape_manual(values = c(19, 17, 15, 18, 4))

# Plots
img(width = 2 * wid, height = wid * 3.5 / 4,
    file = paste0(figloc, 'fig-agg-ts', suf))
print(global.plot + theme_bw())
dev.off()

# Conjoin wedge proportions for global statistics on agg, pure, & applied
apa.dat <- do.call(rbind, list(
    data.frame(Network = 'Aggregate',
               Variety = 'Classical',
               Year = years,
               Statistic = sapply(mathrev.wedges, wedges2cc, 1, 2))
    , data.frame(Network = 'Aggregate',
                 Variety = 'Opsahl',
                 Year = years,
                 Statistic = sapply(mathrev.wedges, wedges2cc, 3, 4))
    , data.frame(Network = 'Aggregate',
                 Variety = 'Exclusive',
                 Year = years,
                 Statistic = sapply(mathrev.wedges, wedges2cc, 5, 6))
    , data.frame(Network = 'Aggregate',
                 Variety = 'BipartiteCorrected',
                 Year = years,
                 Statistic = sapply(mathrev.wedges, wedges2cc, 1, 2) /
                     sapply(mathrev.degree[[1]][-length(mathrev.degree[[1]])],
                            function(lst) {
                                nsw.clust(lst[[1]], lst[[2]])
                            }))
    , data.frame(Network = 'Aggregate',
                 Variety = 'Dynamic',
                 Year = years,
                 Statistic = sapply(mathrev.closed, wedges2tr, 1, 2))
    #, data.frame(Network = 'Aggregate',
    #             Variety = 'Dynamic2',
    #             Year = years,
    #             Statistic = sapply(mathrev.closes[[2]], wedges2tr, 1, 2))
    , data.frame(Network = 'Pure',
                 Variety = 'Classical',
                 Year = years,
                 Statistic = sapply(mathrev.pure.wedges, wedges2cc, 1, 2))
    , data.frame(Network = 'Pure',
                 Variety = 'Opsahl',
                 Year = years,
                 Statistic = sapply(mathrev.pure.wedges, wedges2cc, 3, 4))
    , data.frame(Network = 'Pure',
                 Variety = 'Exclusive',
                 Year = years,
                 Statistic = sapply(mathrev.pure.wedges, wedges2cc, 5, 6))
    , data.frame(Network = 'Pure',
                 Variety = 'BipartiteCorrected',
                 Year = years,
                 Statistic = sapply(mathrev.pure.wedges, wedges2cc, 1, 2) /
                     sapply(mathrev.degree[[2]][-length(mathrev.degree[[2]])],
                            function(lst) {
                                nsw.clust(lst[[1]], lst[[2]])
                            }))
    , data.frame(Network = 'Pure',
                 Variety = 'Dynamic',
                 Year = years,
                 Statistic = sapply(mathrev.pure.closed, wedges2tr, 1, 2))
    #, data.frame(Network = 'Pure',
    #             Variety = 'Dynamic2',
    #             Year = years,
    #             Statistic = sapply(mathrev.pure.closes[[2]], wedges2tr, 1, 2))
    , data.frame(Network = 'Applied',
                 Variety = 'Classical',
                 Year = years,
                 Statistic = sapply(mathrev.applied.wedges, wedges2cc, 1, 2))
    , data.frame(Network = 'Applied',
                 Variety = 'Opsahl',
                 Year = years,
                 Statistic = sapply(mathrev.applied.wedges, wedges2cc, 3, 4))
    , data.frame(Network = 'Applied',
                 Variety = 'Exclusive',
                 Year = years,
                 Statistic = sapply(mathrev.applied.wedges, wedges2cc, 5, 6))
    , data.frame(Network = 'Applied',
                 Variety = 'BipartiteCorrected',
                 Year = years,
                 Statistic = sapply(mathrev.applied.wedges, wedges2cc, 1, 2) /
                     sapply(mathrev.degree[[3]][-length(mathrev.degree[[3]])],
                            function(lst) {
                                nsw.clust(lst[[1]], lst[[2]])
                            }))
    , data.frame(Network = 'Applied',
                 Variety = 'Dynamic',
                 Year = years,
                 Statistic = sapply(mathrev.applied.closed,
                                    wedges2tr, 1, 2))
    #, data.frame(Network = 'Applied',
    #             Variety = 'Dynamic2',
    #             Year = years,
    #             Statistic = sapply(mathrev.applied.closes[[2]],
    #                                wedges2tr, 1, 2))
))

# Restrict to adjacent 3-year intervals
apa.dat <- apa.dat[apa.dat$Year %in% seq(1987, 2008, 3), ]

# Insert empty facet at position 4
apa.dat$Variety.alt <- factor(apa.dat$Variety,
                              levels = c(levels(apa.dat$Variety)[1:3],
                                         '',
                                         levels(apa.dat$Variety)[4:5]))

# Time series plots for global statistics on aggregate, pure, & applied networks
apa.plot <- ggplot(data = apa.dat,
                  aes(x = Year, y = Statistic, shape = Network)) +
    geom_point() + geom_line() +
    facet_wrap(~ Variety.alt, ncol = 3, drop = FALSE, scales = 'free_y') +
    scale_x_continuous('End year') +
    theme_bw() +
    theme(legend.position = c(0, 0), legend.justification = c(0, 0))

# Remove shaded backdrop from this facet
apa.grob <- ggplotGrob(apa.plot)
## which elements to remove
rm.str <- '(panel|strip_t).*4'
## remove empty panels
apa.grob$grobs[grep(rm.str, names(apa.grob$grobs))] <- NULL
## remove them from the layout
apa.grob$layout = apa.grob$layout[-grep(rm.str, (apa.grob$layout)$name), ]

## move bottom axis to above spot
apa.grob$layout[grep('axis_b.*4', apa.grob$layout$name), c('t', 'b')] <-
    apa.grob$layout[grep('axis_b.*4', apa.grob$layout$name), c('t', 'b')] / 1.5
## move left axis to next spot
apa.grob$layout[grep('axis_l.*4', apa.grob$layout$name), c('l', 'r')] <-
    apa.grob$layout[grep('axis_l.*4', apa.grob$layout$name), c('l', 'r')] * 2

# Plots
img(width = 2 * wid, height = wid * 5 / 4,
    file = paste0(figloc, 'fig-apa-ts', suf))
grid.draw(apa.grob)
dev.off()

rm(list = ls())
