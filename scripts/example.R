# List of example networks to use in analysis
library(bitriad)
data(list = data(package = 'bitriad')$results[, 3])

example <- list(DG1 = women_group,
                DG2 = women_clique,
                BB = chicago1960s,
                GWF = minneapolis1970s,
                NMT1 = nmt_meetings,
                NMT2 = nmt_organizations,
                FH = whigs)

save(example, file = 'calc/example.RData')

rm(list = ls())
