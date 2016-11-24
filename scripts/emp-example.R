# Setup
library(bitriad)
library(xtable)
library(ggplot2)
library(reshape2)
#detach("package:network", unload = TRUE) # needs an if-then shell
source('code/triadic-base.R')
source('code/triadic-cent.R')
source('code/bonacich.centrality.R')
source('code/random-bipartite.R')
load('calc/example.RData')
load('calc/example-model.RData')

# Function to evaluate global clustering coefficients
eval.global <- function(bigraph) {
    c(transitivity(actor.projection(bigraph)),
      opsahl.transitivity(bigraph),
      excl.transitivity(bigraph))
}

# Function to evaluate local clustering coefficients
eval.local <- function(bigraph, vids) {
    proj.vids <- which(which(!V(bigraph)$type) %in% vids)
    stopifnot(length(proj.vids) == length(vids))
    el <- cbind(transitivity(actor.projection(bigraph),
                       type = 'local', vids = proj.vids),
          opsahl.transitivity(bigraph, type = 'local', vids = vids),
          excl.transitivity(bigraph, type = 'local', vids = vids))
    if (is.dyn(bigraph)) {
        el <- cbind(el,
                    dyn.transitivity.an(bigraph, type = 'local')[proj.vids])
    }
    el
}

# Identify desired networks
wh.ex <- which(names(example) %in% c('DG1', 'GWF'))

# Compute global statistics
example.gcc <- t(sapply(example[wh.ex], eval.global))
colnames(example.gcc) <- c('Classical', 'Opsahl', 'Exclusive')

# Print global clustering coefficients for each network
print('Global clustering coefficients for small networks')
print(round(example.gcc, 3))

# Compute global and local statistics on DG1 and GWF
example.lcc <- lapply(example[wh.ex], function(bigraph) {
    vids <- local.reps(bigraph)
    mat <- eval.local(bigraph, vids = vids)
    colnames(mat) <-
        c('Classical', 'Opsahl', 'Exclusive', 'Dynamic')[1:ncol(mat)]
    rownames(mat) <- names(vids)
    mat
})

# Compute degree-corrected eigenvector centralities of actors
example.cent <- lapply(example[wh.ex], function(bigraph) {
    # Degree centrality
    deg <- degree(actor.projection(bigraph))
    # Multidegree (2-walk) centrality
    mdeg.cent <- unit.cent(bigraph, 1)
    # 4-walk centrality
    p.cent <- unit.cent(bigraph, 2)
    # Eigenvector centrality
    eigen.cent <- ev.cent(bigraph)
    cent <- data.frame(
        #Degree = deg / sqrt(sum(deg ^ 2)),
        TwoWalk = mdeg.cent / sqrt(sum(mdeg.cent ^ 2)),
        #FourWalk = p.cent / sqrt(sum(p.cent ^ 2)),
        Eigenvector = eigen.cent / sqrt(sum(eigen.cent ^ 2)),
        #DegreeCorrected = eigen.cent / sqrt(sum(eigen.cent ^ 2)) -
        #    mdeg.cent / sqrt(sum(mdeg.cent ^ 2)),
        TwoWalkCorrected = eigen.cent / sqrt(sum(eigen.cent ^ 2)) -
            mdeg.cent / sqrt(sum(mdeg.cent ^ 2))#,
        #FourWalkCorrected = eigen.cent / sqrt(sum(eigen.cent ^ 2)) -
        #    p.cent / sqrt(sum(p.cent ^ 2))
    )
    cent[local.reps(bigraph), ]
})

# Conjoin local clustering coefficients with eigenvector centralities
example.mats <- lapply(1:length(example.lcc), function(i) {
    cbind(example.lcc[[i]], example.cent[[i]])
})
names(example.mats) <- names(example)[wh.ex]

# Write tables of local diagnostics
for (i in 1:length(example.mats)) {
    latex.table(round(example.mats[[i]], 3),
                digits = 3,
                align = paste0('l',
                               paste(rep('r', ncol(example.mats[[i]])),
                                     collapse = '')),
                math.mode = c(),
                file = paste0(tabloc, 'tab-', names(example.mats)[i], '.txt'))
}

# First remove any unshared columns
intersection <- function(lst) {
    if (length(lst) == 1) return(lst[[1]])
    if (length(lst) == 2) return(intersect(lst[[1]], lst[[2]]))
    intersect(lst[[1]], intersection(lst[-1]))
}
incl.col <- intersection(lapply(example.mats, colnames))
# Single data frame for plotting
example.dat <- do.call(rbind, lapply(1:length(example.mats), function(i) {
    data.frame(Network = names(example.mats)[i],
               example.mats[[i]][, incl.col])
}))
example.melt <- melt(example.dat, measure.vars = 2:4) # HOW MANY OF EACH?
example.melt2 <- melt(example.melt, measure.vars = 2:4) # HOW MANY OF EACH?
names(example.melt2)[2:5] <- c('Variety',
                               'ClusteringCoefficient',
                               'Type',
                               'Centrality')

# Plot clustering coefficients versus centralities by network and statistic
cent.plot <- if (cvcc) {
    ggplot(data = example.melt2,
           aes(x = ClusteringCoefficient, y = Centrality))
} else {
    ggplot(data = example.melt2,
           aes(x = Centrality, y = ClusteringCoefficient))
}
cent.plot <- cent.plot +
    geom_point() +
    facet_wrap(Variety + Type ~ Network, ncol = 2, scales = 'free') +
    geom_smooth(method = 'lm', colour = 'black')

# Print plots
img(width = 1.5 * wid, height = 6 * wid * 4 / 4,
    file = paste0(figloc, 'fig-ex-cent-all', suf))
print(cent.plot + theme_bw())
dev.off()

# Which clustering coefficients and centrality scores to include?
wh.clust <- c("Classical", 'Opsahl', 'Exclusive')
wh.cent <- c('TwoWalk', 'TwoWalkCorrected')

# Plot exclusive clustering coefficient versus degree-corrected centrality
cent.subplot <- lapply(levels(example.melt2$Network), function(ntwk) {
    subplot <- if (cvcc) {
        ggplot(
            data = example.melt2[example.melt2$Variety %in% wh.clust &
                                     example.melt2$Type %in% wh.cent &
                                     example.melt2$Network == ntwk, ],
            aes(x = ClusteringCoefficient, y = Centrality)
        )
    } else {
        ggplot(
            data = example.melt2[example.melt2$Variety %in% wh.clust &
                                     example.melt2$Type %in% wh.cent &
                                     example.melt2$Network == ntwk, ],
            aes(x = Centrality, y = ClusteringCoefficient)
        )
    }
    subplot <- subplot +
        geom_point() +
        facet_wrap(~ Type + Variety, scales = 'free', ncol = 3) +
        geom_smooth(method = 'lm', colour = 'black')
    subplot <- if (cvcc) {
        subplot + xlab('Clustering coefficient')
    } else {
        subplot + ylab('Clustering coefficient')
    }
    subplot
})
names(cent.subplot) <- levels(example.melt2$Network)

# Print subplot
for (ntwk in levels(example.melt2$Network)) {
    img(width = 2 * wid, height = wid * 4 / 4,
        file = paste0(figloc, 'fig-ex-cent-', ntwk, suf))
    print(cent.subplot[[ntwk]] + theme_bw())
    dev.off()
}

# Linear models of clustering coefficients from degree-corrected centralities
cent.model <- lapply(unique(example.melt2$Network), function(ntwk) {
    models <- lapply(wh.clust, function(vr) {
        submodels <- lapply(wh.cent, function(ty) {
            dat <- example.melt2[example.melt2$Network == ntwk &
                                    example.melt2$Variety == vr &
                                    example.melt2$Type == ty, ]
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
names(cent.model) <- unique(example.melt2$Network)

# Summaries of models
print('Clustering coefficient-centrality regression summaries')
for (ntwk in unique(example.melt$Network)) {
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

# Set number of digits
ndig <- 3



# Rename model columns
for (i in 1:length(example.model)) {
    for (j in 1:length(example.model[[i]])) {
        names(example.model[[i]][[j]])[1:3] <-
            c('Classical',
              'Opsahl',
              'Exclusive')
    }
}

# Subset of networks to include in global null model plot
example.incl <- c('DG1', 'GWF')
index.incl <- which(names(example.model) %in% example.incl)
global.lst <- lapply(example.model[index.incl], function(lst) lst$global)

# Global null model plot
img(width = 2 * wid, height = length(index.incl) * wid * 1 / 3,
    file = paste0(figloc, 'fig-ex-null', suf))
hist.grid(lst = global.lst, col = 'darkgray', lty = 'dashed')
dev.off()

# Local null model plots
for (i in 1:length(example.incl)) {
    
    # Identify the model results for the desired network's structural reps
    reps <- local.reps(example[[example.incl[[i]]]])
    local.lst <- example.model[[example.incl[[i]]]][1 + reps]
    
    # Grid of histograms!
    img(width = 2 * wid, height = 12,
        file = paste0(figloc, 'fig-', example.incl[i], '-null', suf))
    hist.grid(lst = local.lst, col = 'darkgray', lty = 'dotted')
    dev.off()
    
}

rm(list = ls())
