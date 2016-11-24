# Setup
library(Matrix)
library(bitriad)
library(ggplot2)
library(reshape2)
#detach("package:network", unload = TRUE) # needs an if-then shell
source('code/triadic-base.R')
source('code/triadic-cent.R')
#source('code/random-bipartite.R')

load('calc/mathrev-wedges.RData')
load("calc/mathrev-center.RData")

# Convert wedge data frames to local clustering coefficient data frames
mathrev.lcc <- lapply(inds, function(i) {
    data.frame(
        Classical = mathrev.wedges[[i]][, 2] / mathrev.wedges[[i]][, 1],
        Opsahl = mathrev.wedges[[i]][, 4] / mathrev.wedges[[i]][, 3],
        Exclusive = mathrev.wedges[[i]][, 6] / mathrev.wedges[[i]][, 5]
    )
})

# Convert pure wedge data frames to local clustering coefficient data frames
mathrev.pure.lcc <- lapply(inds, function(i) {
    data.frame(
        Classical = mathrev.pure.wedges[[i]][, 2] /
            mathrev.pure.wedges[[i]][, 1],
        Opsahl = mathrev.pure.wedges[[i]][, 4] /
            mathrev.pure.wedges[[i]][, 3],
        Exclusive = mathrev.pure.wedges[[i]][, 6] /
            mathrev.pure.wedges[[i]][, 5]
    )
})

# Convert applied wedge data frames to local clustering coefficient data frames
mathrev.applied.lcc <- lapply(inds, function(i) {
    data.frame(
        Classical = mathrev.applied.wedges[[i]][, 2] /
            mathrev.applied.wedges[[i]][, 1],
        Opsahl = mathrev.applied.wedges[[i]][, 4] /
            mathrev.applied.wedges[[i]][, 3],
        Exclusive = mathrev.applied.wedges[[i]][, 6] /
            mathrev.applied.wedges[[i]][, 5]
    )
})

# Conjoin local clustering coefficients with eigenvector centralities
mathrev.mats <- lapply(1:length(mathrev.cent), function(i) {
    cbind(mathrev.lcc[[i]], mathrev.cent[[i]])
})
names(mathrev.mats) <- yrs

# Single data frame for plotting
mathrev.dat <- do.call(rbind, lapply(1:length(mathrev.mats), function(i) {
    data.frame(Network = paste0("MR (", names(mathrev.mats)[i], ")"),
               mathrev.mats[[i]])
}))
mathrev.melt <- melt(mathrev.dat, measure.vars = 2:4) # HOW MANY OF EACH?
mathrev.melt2 <- melt(mathrev.melt, measure.vars = 2:4) # HOW MANY OF EACH?
names(mathrev.melt2)[2:5] <- c('Variety',
                               'ClusteringCoefficient',
                               'Type',
                               'Centrality')

# Plot clustering coefficients versus centralities by network and statistic
cent.plot <- if (cvcc) {
    ggplot(data = mathrev.melt2,
           aes(x = ClusteringCoefficient, y = Centrality))
} else {
    ggplot(data = mathrev.melt2,
           aes(x = Centrality, y = ClusteringCoefficient))
}
cent.plot <- cent.plot +
    geom_point() +
    facet_wrap(Variety + Type ~ Network,
               ncol = length(levels(mathrev.melt2$Network)),
               scales = 'free') +
    geom_smooth(method = 'lm', colour = 'black')

# Print plots
img(width = 2 * wid, height = 6 * wid * 4 / 4,
    file = paste0(figloc, 'fig-mr-cent-all', suf))
print(cent.plot + theme_bw())
dev.off()

# Which clustering coefficients and centrality scores to include?
wh.clust <- c("Classical", 'Opsahl', 'Exclusive')
wh.cent <- c('TwoWalk', 'TwoWalkCorrected')

# Plot exclusive clustering coefficient versus degree-corrected centrality
cent.subplot <- lapply(levels(mathrev.melt2$Network), function(ntwk) {
    subplot <- if (cvcc) {
        ggplot(
            data = mathrev.melt2[mathrev.melt2$Variety %in% wh.clust &
                                     mathrev.melt2$Type %in% wh.cent &
                                     mathrev.melt2$Network == ntwk, ],
            aes(x = ClusteringCoefficient, y = Centrality)
        )
    } else {
        ggplot(
            data = mathrev.melt2[mathrev.melt2$Variety %in% wh.clust &
                                     mathrev.melt2$Type %in% wh.cent &
                                     mathrev.melt2$Network == ntwk, ],
            aes(x = Centrality, y = ClusteringCoefficient)
        )
    }
    subplot <- subplot +
        geom_point(colour = rgb(0, 0, 0, .1)) +
        facet_wrap(~ Type + Variety, scales = 'free',
                   ncol = length(levels(mathrev.melt2$Variety))) +
        geom_smooth(method = 'lm', colour = 'black')
    subplot <- if (cvcc) {
        subplot + xlab('Clustering coefficient')
    } else {
        subplot + ylab('Clustering coefficient')
    }
    subplot
})
names(cent.subplot) <- levels(mathrev.melt2$Network)

# Print subplot
for (ntwk in levels(mathrev.melt2$Network)) {
    img(width = 2 * wid, height = wid * 4 / 4,
        file = paste0(figloc, 'fig-mr-cent-', ntwk, suf))
    print(cent.subplot[[ntwk]] + theme_bw())
    dev.off()
}

# Linear models of clustering coefficients from degree-corrected centralities
cent.model <- lapply(unique(mathrev.melt2$Network), function(ntwk) {
    models <- lapply(wh.clust, function(vr) {
        submodels <- lapply(wh.cent, function(ty) {
            dat <- mathrev.melt2[mathrev.melt2$Network == ntwk &
                                     mathrev.melt2$Variety == vr &
                                     mathrev.melt2$Type == ty, ]
            if (cvcc) {
                lm(data = dat,
                   formula = Centrality ~ ClusteringCoefficient)
            } else {
                lm(data = dat,
                   formula = ClusteringCoefficient ~ Centrality)
            }
            
        })
        names(submodels) <- wh.cent
        submodels
    })
    names(models) <- wh.clust
    models
})
names(cent.model) <- unique(mathrev.melt2$Network)

# Summaries of models
print('Clustering coefficient-centrality regression summaries')
for (ntwk in unique(mathrev.melt$Network)) {
    for (vr in wh.clust) {
        for (ty in wh.cent) {
            print(paste(vr,
                        if (cvcc) 'centrality' else 'clustering',
                        'regressed on', ty,
                        if (cvcc) 'clustering' else 'cetrality',
                        'in', ntwk))
            print(summary(cent.model[[ntwk]][[vr]][[ty]]))
        }
    }
}

rm(list = ls())
