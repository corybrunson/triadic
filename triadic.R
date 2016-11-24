# Install necessary packages
ipak <- function(pkg) {
    new.pkg <- pkg[!(pkg %in% installed.packages()[, 'Package'])]
    if(length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
}
ipak(c('devtools',
       'igraph', 'networksis',
       'ggplot2', 'xtable', 'grid',
       'reshape2', 'dplyr'))
for(pkg in c('bitriad')) {
    if(!(pkg %in% installed.packages()[, 'Package']))
        install_github(paste0('corybrunson/', pkg))
}

# Create directories
if(grepl('/triadic$', getwd())) {
    dirs <- c('calc', 'fig', 'tab')
    for(d in dirs) if(!file.exists(d)) dir.create(d)
} else stop('Go to (or create) directory "triadic"')

rm(list = ls())

# Clean data sets
source('scripts/example.R') # single data frame
source('scripts/mathrev.R') # separate graph objects

# Run intensive calculations from which tables and figures can be derived
source('scripts/example-census.R')
source('scripts/example-models.R')
source('scripts/mathrev-degree.R')
source('scripts/mathrev-census.R')
source('scripts/mathrev-wedges.R')
source('scripts/mathrev-closed.R')
source('scripts/mathrev-closes.R')
source('scripts/mathrev-center.R')
source('scripts/mathrev-char2s.R')

# Illustrative examples
source('scripts/fig-kite.R')
source('scripts/fig-dg2.R')
source('scripts/tab-dg2.R')
source('scripts/fig-triads.R')
source('scripts/fig-schedule.R')

# Empirical tests
source('scripts/emp-censuses.R')
source('scripts/emp-example.R')
source('scripts/emp-test.R')
source('scripts/emp-apa.R')
source('scripts/emp-stc.R')
source('scripts/emp-dep.R')
source('scripts/emp-cent.R')
