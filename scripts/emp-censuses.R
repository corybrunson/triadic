# Setup
library(bitriad)
library(xtable)
source('code/triadic-base.R')
load('calc/example-census.RData')
load('calc/mathrev-census.RData')

# Which MR intervals to include
int.incl <- c(1, 3)

# DG1, GWF, and two MR structural censuses
ftc2stc <- function(ftc) project_census(ftc, scheme = "full")$binary
stc.list <- c(
    lapply(example.census[c('DG1', 'GWF')], ftc2stc)
    , lapply(mathrev.census[int.incl], ftc2stc)
)
names(stc.list) <- c(
    c('DG1', 'GWF')
    , paste0('MR (',
             paste0(yrs[int.incl] - 2, '-', substr(yrs[int.incl], 4, 4)),
             ')')
)

# Write tables together
latex.tables(list = stc.list, digits = 0,
             align.row = 'c|', align.each = 'rr|', align.headers = 'c|',
             file = paste0(tabloc, 'tab-ex-census.txt'),
             math.mode = 'entries')

rm(list = ls())
