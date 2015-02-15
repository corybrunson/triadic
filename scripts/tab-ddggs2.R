# Table: Full census and triadic closure measures for DDGGS2
library(bitriad)
library(xtable)
source('code/triadic-base.R')
load('calc/example.RData')
load('calc/example-census.RData')

# Full triad census
ddggs2.ftc <- example.census$DDGGS2
rownames(ddggs2.ftc) <- paste0('(',
                          apply(sapply(0:(nrow(ddggs2.ftc) - 1),
                                       index.partition),
                                2, paste, collapse = ','),
                          ')')
colnames(ddggs2.ftc) <- as.character(0:(ncol(ddggs2.ftc) - 1))
# Write to a LaTeX table
latex.table(ddggs2.ftc, digits = 0,
            align = paste0('c|',
                           paste(rep('r', ncol(ddggs2.ftc)), collapse = ''),
                           '|'),
            file = paste0(tabloc, 'tab-ddggs2-census.txt'), math.mode = 'all')
# Write to a matrix with LaTeX math formatting
ddggs2.ftc.mat <- cbind(
    rownames(ddggs2.ftc),
    ddggs2.ftc
)
write(paste0(' & ',
             paste(colnames(ddggs2.ftc), collapse = ' & '),
             '\\\\\\hline\n'),
      file = paste0(tabloc, 'tab-ddggs2-census.txt'))
write.table(ddggs2.ftc.mat,
            quote = FALSE, sep = ' & ', na = '', eol = ' \\\\\n',
            row.names = FALSE, col.names = FALSE,
            file = paste0(tabloc, 'tab-ddggs2-census.txt'), append = TRUE)

# Time stamps
ddggs2 <- example$DDGGS2
V(ddggs2)$time[V(ddggs2)$type] <- 1:5

# Measures of triadic closure
ddggs2.tc <- rbind(
    c(
        transitivity(actor.projection(ddggs2), type = 'global')
        , opsahl.transitivity(ddggs2, type = 'global')
        , excl.transitivity(ddggs2, type = 'global')
        #, dyn.triadic.closure(ddggs2, type = 'global')
    ),
    do.call(cbind, list(
        transitivity(actor.projection(ddggs2), type = 'local')
        , opsahl.transitivity(ddggs2, type = 'local')
        , excl.transitivity(ddggs2, type = 'local')
        #, dyn.triadic.closure(ddggs2, type = 'local')
    ))
)

# Row and column names
rownames(ddggs2.tc) <- c('DDGGS2', V(ddggs2)$name[!V(ddggs2)$type])
colnames(ddggs2.tc) <- c(
    'Classical'
    , 'Opsahl'
    , 'Exclusive'
    #, 'Dynamic'
)

# Set number of digits and replace undefined entries with empty strings
ndig <- 3
wh.nan <- which(is.nan(ddggs2.tc))
ddggs2.tc <- format(round(ddggs2.tc, ndig), nsmall = ndig)
ddggs2.tc[wh.nan] <- ''

# Write table
latex.table(ddggs2.tc, digits = 0,
            align = paste0('c',
                           paste(rep('r', ncol(ddggs2.tc)), collapse = ''),
                           ''),
            file = paste0(tabloc, 'tab-ddggs2-tc.txt'), math.mode = 'entries')

rm(list = ls())
