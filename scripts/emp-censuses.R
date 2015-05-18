# Setup
library(bitriad)
library(xtable)
source('code/triadic-base.R')
load('calc/example-census.RData')
load('calc/mathrev-census.RData')

# Which MR intervals to include
int.incl <- c(1, 3)

# DDGGS2, DDGGS1, GWF, and two MR structural censuses
ftc2stc <- function(ftc) project.census(ftc, scheme = "full")
stc.list <- c(
    lapply(example.census[c('DDGGS2', 'DDGGS1', 'GWF')], ftc2stc)
    , lapply(mathrev.census[int.incl], ftc2stc)
)
names(stc.list) <- c(
    c('DDGGS2', 'DDGGS1', 'GWF')
    , paste0('MR (',
             paste0(yrs[int.incl] - 2, '-', substr(yrs[int.incl], 4, 4)),
             ')')
)

# Row and column names
for(i in 1:length(stc.list)) {
    rownames(stc.list[[i]]) <- 0:3
    colnames(stc.list[[i]]) <- 0:1
}

# Write tables together
latex.tables(list = stc.list, digits = 0,
             align.row = 'c|', align.each = 'rr|', align.headers = 'c|',
             file = paste0(tabloc, 'tab-ex-census.txt'),
             math.mode = 'entries')

rm(list = ls())
