# List of example networks to use in analysis
library(bitriad)
data(list = data(package = 'bitriad')$results[, 3])

example <- list(DDGGS1 = ddggs.group,
                DDGGS2 = ddggs.clique,
                BB = barnes.burkett.corporate,
                GWF = galaskiewicz.ceos,
                NMT1 = noordin.top.meetings,
                NMT2 = noordin.top.organizations,
                FH = fischer.han.whigs)

save(example, file = 'calc/example.RData')

rm(list = ls())
