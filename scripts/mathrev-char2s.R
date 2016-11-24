# Figure: Overlapping time series of global diagnostics on the MR network
library(igraph)
library(bitriad)
source('code/mathrev2igraph.R')
source('code/triadic-spec.R')
load('calc/mathrev.RData')

# List of prefixes
char2s <- sort(unique(substr(mathrev$pclass, 1, 2)))
char2s <- char2s[grep('[0-9]{2}', char2s)]

# Degree sequences
char2s.degree <- lapply(years[-(1:2)], function(yr) {
  dat <- mathrev[mathrev$year %in% (yr - dur + 1):yr, ]
  lst <- lapply(char2s, function(s) {
    bigraph <- as_an(paper.author.graph(
      dat[substr(dat$pclass, 1, 2) == s, ]
    ))
    if(vcount(bigraph) == 0) return(list(actor = c(), event = c()))
    deg <- degree(bigraph)
    list(
      actor = tabulate(deg[which(!V(bigraph)$type)]),
      event = tabulate(deg[which(V(bigraph)$type)])
    )
  })
})

# Compute triad censuses for each subdiscipline
# over three evenly-spaced intervals of equal duration
char2s.census <- lapply(years[-(1:2)], function(yr) {
  dat <- mathrev[mathrev$year %in% (yr - dur + 1):yr, ]
  lst <- lapply(char2s, function(s) {
    triad_census_an(as_an(paper.author.graph(
      dat[substr(dat$pclass, 1, 2) == s, ])))
  })
  names(lst) <- char2s
  lst
})

# Compute global statistics for each subdiscipline
char2s.global <- lapply(char2s.census, function(lst) {
  do.call(rbind, lapply(lst, function(census) {
    c("C" = transitivity_from_census(census, "classical", scheme = "full"),
      "C.opsahl" = transitivity_from_census(census, "opsahl", scheme = "full"),
      "C.excl" = transitivity_from_census(census, "exclusive", scheme = "full"),
      "C.injstr" = transitivity_from_census(census, "injstr", scheme = "full"),
      "C.injact" = transitivity_from_census(census, "injact", scheme = "full"),
      "T" = transitivity_from_census(census, "classical", scheme = "full"),
      "T.excl" = transitivity_from_census(census, "exclusive", scheme = "full"),
      "T.injact" = transitivity_from_census(census, "injact", scheme = "full")
    )
  }))
})

# Compute wedge closure (both types) for each subdiscipline
char2s.closes <- lapply(years[-(1:2)], function(yr) {
  dat <- mathrev[mathrev$year %in% (yr - dur + 1):yr, ]
  do.call(rbind, lapply(char2s, function(s) {
    sdat <- dat[substr(dat$pclass, 1, 2) == s, ]
    return(c(DTC1 = wedge.closure(sdat, (yr - dur + 1), (yr - dur + 2):yr),
             DTC2 = wedge.closure(sdat, (yr - dur + 1):(yr - dur + 2), yr)))
  }))
})

# Compute wedge closure (each type) for each subdiscipline
char2s.closed <- lapply(years[-(1:2)], function(yr) {
  dat <- mathrev[mathrev$year %in% (yr - dur + 1):yr, ]
  sapply(char2s, function(s) {
    sdat <- dat[substr(dat$pclass, 1, 2) == s, ]
    dynamic_transitivity_an(as_an(paper.author.graph(sdat)),
                            memory = Inf,
                            type = 'global')
  })
})

# Save!
save(char2s.degree, char2s.census, char2s.global, char2s.closes, char2s.closed,
     file = 'calc/mathrev-char2s.RData')

rm(list = ls())
