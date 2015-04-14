# List of example networks to use in analysis
library(bitriad)
data(list = data(package = 'bitriad')$results[, 3])

example <- list(DDGGS1 = davis.group,
                DDGGS2 = davis.clique,
                BB = chicago1960s,
                GWF = minneapolis1970s,
                NMT1 = nmt.meetings,
                NMT2 = nmt.organizations,
                FH = whigs)

save(example, file = 'calc/example.RData')

rm(list = ls())
