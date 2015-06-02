# Table: Full census and triadic closure measures for DG2
library(bitriad)
library(xtable)
source('code/triadic-base.R')
load('calc/example.RData')
load('calc/example-census.RData')

# Full triad census
dg2.ftc <- example.census$DG2
rownames(dg2.ftc) <- paste0('(',
                          apply(sapply(0:(nrow(dg2.ftc) - 1),
                                       indexPartition),
                                2, paste, collapse = ','),
                          ')')
colnames(dg2.ftc) <- as.character(0:(ncol(dg2.ftc) - 1))
# Write to a LaTeX table
latex.table(dg2.ftc, digits = 0,
            align = paste0('c|',
                           paste(rep('r', ncol(dg2.ftc)), collapse = ''),
                           '|'),
            file = paste0(tabloc, 'tab-dg2-census.txt'), math.mode = 'all')
# Write to a matrix with LaTeX math formatting
dg2.ftc.mat <- cbind(
    rownames(dg2.ftc),
    dg2.ftc
)
write(paste0(' & ',
             paste(colnames(dg2.ftc), collapse = ' & '),
             '\\\\\\hline\n'),
      file = paste0(tabloc, 'tab-dg2-census.txt'))
write.table(dg2.ftc.mat,
            quote = FALSE, sep = ' & ', na = '', eol = ' \\\\\n',
            row.names = FALSE, col.names = FALSE,
            file = paste0(tabloc, 'tab-dg2-census.txt'), append = TRUE)

# Time stamps
dg2 <- example$DG2
V(dg2)$time[V(dg2)$type] <- 1:5

# Measures of triadic closure
dg2.tc <- rbind(
    c(
        transitivity(actor.projection(dg2), type = 'global')
        , opsahl.transitivity(dg2, type = 'global')
        , excl.transitivity(dg2, type = 'global')
        #, dyn.transitivity(dg2, type = 'global')
    ),
    do.call(cbind, list(
        transitivity(actor.projection(dg2), type = 'local')
        , opsahl.transitivity(dg2, type = 'local')
        , excl.transitivity(dg2, type = 'local')
        #, dyn.transitivity(dg2, type = 'local')
    ))
)

# Row and column names
rownames(dg2.tc) <- c('DG2', V(dg2)$name[!V(dg2)$type])
colnames(dg2.tc) <- c(
    'Classical'
    , 'Opsahl'
    , 'Exclusive'
    #, 'Dynamic'
)

# Set number of digits and replace undefined entries with empty strings
ndig <- 3
wh.nan <- which(is.nan(dg2.tc))
dg2.tc <- format(round(dg2.tc, ndig), nsmall = ndig)
dg2.tc[wh.nan] <- ''
dg2.tc <- t(dg2.tc)

# Write table
latex.table(dg2.tc, digits = 0,
            align = paste0('c',
                           paste(rep('r', ncol(dg2.tc)), collapse = ''),
                           ''),
            file = paste0(tabloc, 'tab-dg2-tc.txt'), math.mode = 'entries')

rm(list = ls())
