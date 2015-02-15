# Setup
library(bitriad)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(schoolmath)
source('code/triadic-base.R')
source('code/random.bipartite.R')
load('calc/example.RData')
load('calc/example-model.RData')
load('calc/mathrev-char2s.RData')

LCM <- function(vec) {
    if(length(vec) == 2) scm(vec[1], vec[2]) else scm(vec[1], LCM(vec[-1]))
}

# Clustering coefficients to include in the analysis
stat.incl <- c('Classical', 'Opsahl', 'Exclusive')
ex.incl <- c('DDGGS1', 'GWF')

# Conjoin static global statistics with dynamic triadic closure
char2s.mats <- lapply(1:length(char2s.global), function(i) {
    cbind(
        char2s.global[[i]][, 1:3]
        , DTC = char2s.closed[[i]]
        #, char2s.closes[[i]]
    )
})

# Professional names
for(i in 1:length(char2s.mats)) colnames(char2s.mats[[i]]) <- c(
    stat.incl, 'Dynamic'
)

# Restrict global statistics to adjacent 3-year intervals
wh.inds <- seq(1, length(char2s.mats), 3)

# Identify subnetworks on which all statistics are defined
wh.def <- Reduce(intersect, lapply(char2s.mats[wh.inds], function(ccmat) {
    which(!is.na(rowSums(ccmat)) &
              apply(ccmat, 1, function(x) length(which(x == 0)) <= 1))
}))

# Define sublist of submatrices for analysis
mr.mats <- lapply(char2s.mats[wh.inds], function(mat) mat[wh.def, ])

# Add degree-corrected (classical) clustering coefficient
for(i in 1:length(wh.inds)) {
    mr.mats[[i]] <- cbind(
        mr.mats[[i]],
        DegreeCorrected = mr.mats[[i]][, 'Classical'] /
            sapply(wh.def, function(j) {
                nsw.clust(char2s.degree[[wh.inds[i]]][[j]][[1]],
                          char2s.degree[[wh.inds[i]]][[j]][[2]])
            })
    )
}

# Construct single matrix with years
mr.mat <- do.call(rbind, mr.mats)

# Local evaluations
ex.mats <- lapply(ex.incl, function(name) {
    sapply(c(watts.strogatz.transitivity,
             opsahl.transitivity,
             excl.transitivity), function(tr) {
                 tr(example[[name]], type = 'local')
             })[local.reps(example[[name]]), ]
})
for(i in 1:length(ex.mats)) colnames(ex.mats[[i]]) <- stat.incl
names(ex.mats) <- ex.incl

# Add degree-corrected (classical) clustering coefficient
for(name in ex.incl) {
    ex.mats[[name]] <- cbind(
        ex.mats[[name]],
        DegreeCorrected = ex.mats[[name]][, 'Classical'] /
            sapply(example.model[[name]][-1], function(mat) {
                sum(mat[-1, 'C'] * mat[-1, 'W']) / sum(mat[-1, 'W'])
            })
    )
}

# Add dynamic triadic closure for dynamic networks
for(name in ex.incl) if(is.dynamic(example[[name]])) {
    ex.mats[[name]] <- cbind(
        ex.mats[[name]],
        Dynamic = dyn.triadic.closure(
            example[[name]], type = 'local'
        )[local.reps(example[[name]])]
    )
}

# Test names
test.names <- c(ex.incl, 'MR')

# STABILITY
# Pooled residual standard deviation
# for values of each C on MRs at adjacent intervals
mr.stability.pool <- sapply(stat.incl, function(name) {
    x <- mr.mat[1:(nrow(mr.mat) - nrow(mr.mats[[length(mr.mats)]])), name]
    y <- mr.mat[(nrow(mr.mats[[1]]) + 1):nrow(mr.mat), name]
    dat <- data.frame(
        value = c(x, y),
        number = rep(as.character(1:nrow(mr.mats[[1]])), times = 2)
    )
    anal.var <- summary(aov(data = dat, value ~ number))
    anal.var[[1]][1, 2] / sum(anal.var[[1]][, 2])
})
# BEGIN WRONG STABILITY STATS
if(FALSE) {
    mr.stability.pool <- sapply(stat.incl, function(name) {
        x <- mr.mat[1:(nrow(mr.mat) - nrow(mr.mats[[length(mr.mats)]])), name]
        y <- mr.mat[(nrow(mr.mats[[1]]) + 1):nrow(mr.mat), name]
        yhat <- predict(lm(y ~ x))
        if(yhat[which(x == min(x))[1]] > yhat[which(x == max(x))][1]) {
            warning('Correlation is negative')
        }
        sum((yhat - mean(y)) ^ 2) / sum((y - mean(y)) ^ 2)
    })
}
# END WRONG STABILITY STATS
stability.mat <- rbind(
    matrix(NA, nr = 2, nc = length(mr.stability.pool)),
    mr.stability.pool
)
rownames(stability.mat) <- test.names
# Mean-difference plots
mr.stability.dat <- data.frame(
    Statistic = rep(colnames(mr.mat)[1:3],
                    each = nrow(mr.mat) - nrow(mr.mats[[1]])),
    Present = as.vector(mr.mat[(nrow(mr.mats[[1]]) + 1):nrow(mr.mat), 1:3]),
    Next = as.vector(mr.mat[1:(nrow(mr.mat) - nrow(mr.mats[[1]])), 1:3])
)
mr.stability.dat$Statistic <- factor(mr.stability.dat$Statistic,
                                     colnames(mr.mat)[1:3])
mr.stability.plot <- ggplot(data = mr.stability.dat,
                            aes(x = (Present + Next) / 2,
                                y = Next - Present)) +
    geom_point() +
    geom_abline(slope = 0, intercept = 0) +
    expand_limits(y = 0) +
    facet_wrap(~ Statistic, ncol = 3, scales = 'free') +
    theme_bw()
img(width = 2 * wid, height = 3 / 4 * wid,
    file = paste0(figloc, 'fig-stability', suf))
plot(mr.stability.plot)
dev.off()

# CONCURRENT VALIDITY
# (Pooled) coefficient of determination
# between values of pairs of a C with a stat it is hypothesized to predict
ex.validity.mat <- t(sapply(ex.mats, function(mat) c(
    Classical = NA,
    Opsahl = cor(mat[, 'Opsahl'], mat[, 'DegreeCorrected'], method = 'pearson'),
    Exclusive = if('Dynamic' %in% colnames(mat)) {
        x <- mat[, 'Exclusive']
        y <- mat[, 'Dynamic']
        yhat <- predict(lm(y ~ x))
        if(yhat[which(x == min(x))[1]] > yhat[which(x == max(x))[1]]) {
            warning('Correlation is negative')
        }
        sum((yhat - mean(y)) ^ 2) / sum((y - mean(y)) ^ 2)
    } else NA
)))
mr.validity.pool <- c(
    Classical = NA,
    Opsahl = {
        x <- mr.mat[, 'Opsahl']
        y <- mr.mat[, 'DegreeCorrected']
        yhat <- predict(lm(y ~ x))
        if(yhat[which(x == min(x))[1]] > yhat[which(x == max(x))[1]]) {
            warning('Correlation is negative')
        }
        sum((yhat - mean(y)) ^ 2) / sum((y - mean(y)) ^ 2)
    },
    Exclusive = {
        x <- mr.mat[, 'Exclusive']
        y <- mr.mat[, 'Dynamic']
        yhat <- predict(lm(y ~ x))
        if(yhat[which(x == min(x))[1]] > yhat[which(x == max(x))[1]]) {
            warning('Correlation is negative')
        }
        sum((yhat - mean(y)) ^ 2) / sum((y - mean(y)) ^ 2)
    }
)
validity.mat <- rbind(
    ex.validity.mat,
    mr.validity.pool
)
rownames(validity.mat) <- test.names
# Residual plots
mr.validity.lm <- list(
    lm(mr.mat[, 'DegreeCorrected'] ~ mr.mat[, 'Opsahl']),
    lm(mr.mat[, 'Dynamic'] ~ mr.mat[, 'Exclusive'])
)
mr.validity.dat <- data.frame(
    Statistic = rep(c('Opsahl', 'Exclusive'), each = nrow(mr.mat)),
    Value = as.vector(mr.mat[, c('Opsahl', 'Exclusive')]),
    Residual = c(mr.validity.lm[[1]]$residuals, mr.validity.lm[[2]]$residuals)
)
mr.validity.dat$Statistic <- factor(mr.validity.dat$Statistic,
                                     colnames(mr.mat)[2:3])
mr.validity.plot <- ggplot(data = mr.validity.dat,
                            aes(x = Value, y = Residual)) +
    geom_point() +
    geom_abline(slope = 0, intercept = 0) +
    expand_limits(y = 0) +
    facet_wrap(~ Statistic, ncol = 3, scales = 'free') +
    theme_bw()
img(width = 2 / 3 * 2 * wid, height = 3 / 4 * wid,
    file = paste0(figloc, 'fig-validity', suf))
plot(mr.validity.plot)
dev.off()

# DISTINGUISHABILITY
# (Pooled) coefficient of "indetermination" between pairs of C
ex.distinguishability.mat <- t(sapply(ex.mats, function(mat) {
    sapply(1:(3 - 1), function(i) sapply(1:3, function(j) {
        if(i >= j) NA else {
            x <- mat[, i]
            y <- mat[, j]
            yhat <- predict(lm(y ~ x))
            sum((y - yhat) ^ 2) / sum((y - mean(y)) ^ 2)
        }
    }))
}))
colnames(ex.distinguishability.mat) <- rep(stat.incl, times = 2)
mr.distinguishability.pool <- as.vector(sapply(1:(3 - 1), function(i) {
    sapply(1:3, function(j) {
        if(i >= j) NA else {
            x <- mr.mat[, i]
            y <- mr.mat[, j]
            yhat <- predict(lm(y ~ x))
            sum((y - yhat) ^ 2) / sum((y - mean(y)) ^ 2)
        }
    })
}))
names(mr.distinguishability.pool) <- rep(stat.incl, times = 2)
distinguishability.mat <- rbind(
    ex.distinguishability.mat,
    mr.distinguishability.pool
)
rownames(distinguishability.mat) <- test.names
# Residual plots
mr.distinguishability.lm <- list(
    lm(mr.mat[, 'Opsahl'] ~ mr.mat[, 'Classical']),
    lm(mr.mat[, 'Exclusive'] ~ mr.mat[, 'Classical']),
    lm(mr.mat[, 'Exclusive'] ~ mr.mat[, 'Opsahl'])
)
mr.distinguishability.dat <- data.frame(
    Comparison = rep(c('Classical-Opsahl',
                       'Classical-Exclusive',
                       'Opsahl-Exclusive'), each = nrow(mr.mat)),
    Value = as.vector(mr.mat[, c('Classical', 'Classical', 'Opsahl')]),
    Residual = unlist(lapply(mr.distinguishability.lm,
                             function(fit) fit$residuals))
)
mr.distinguishability.dat$Comparison <-
    factor(mr.distinguishability.dat$Comparison,
           c('Classical-Opsahl',
             'Classical-Exclusive',
             'Opsahl-Exclusive'))
mr.distinguishability.plot <- ggplot(data = mr.distinguishability.dat,
                           aes(x = Value, y = Residual)) +
    geom_point() +
    geom_abline(slope = 0, intercept = 0) +
    expand_limits(y = 0) +
    facet_wrap(~ Comparison, ncol = 3, scales = 'free') +
    theme_bw()
img(width = 2 * wid, height = 3 / 4 * wid,
    file = paste0(figloc, 'fig-distinguishability', suf))
plot(mr.distinguishability.plot)
dev.off()

# DISCRIMINABILITY
# (Pooled) standardized variance of values of each C
ex.discriminability.mat <- t(sapply(ex.mats, function(mat) {
    apply(mat[, 1:3], 2, function(x) {
        (sum((x - mean(x)) ^ 2) / (length(x) - 1)) / .25
    })
}))
mr.discriminability.pool <- sapply(stat.incl, function(name) {
    x <- mr.mat[, name]
    (sum((x - mean(x)) ^ 2) / (length(x) - 1)) / .25
})
discriminability.mat <- rbind(
    ex.discriminability.mat,
    mr.discriminability.pool
)
rownames(discriminability.mat) <- test.names
# Histograms
mr.discriminability.dat <- data.frame(
    Statistic = rep(c('Classical', 'Opsahl', 'Exclusive'), each = nrow(mr.mat)),
    Value = as.vector(mr.mat[, c('Classical', 'Opsahl', 'Exclusive')])
)
mr.discriminability.dat$Statistic <-
    factor(mr.discriminability.dat$Statistic,
           c('Classical', 'Opsahl', 'Exclusive'))
mr.discriminability.plot <- ggplot(data = mr.discriminability.dat,
                                     aes(x = Value)) +
    geom_bar() + stat_bin(binwidth = .02) +
    xlim(0, 1) +
    ylab('Count') +
    facet_wrap(~ Statistic, ncol = 3, scales = 'free') +
    theme_bw()
img(width = 2 * wid, height = 2 / 4 * wid,
    file = paste0(figloc, 'fig-discriminability', suf))
plot(mr.discriminability.plot)
dev.off()

# Bind into a single matrix
test.mat <- rbind(
    as.vector(stability.mat)
    , as.vector(validity.mat)
    , matrix(as.vector(distinguishability.mat), nrow = 2, byrow = TRUE)
    #, rep(NA, length(stat.incl) * length(test.names))
    , as.vector(discriminability.mat)
)
rownames(test.mat) <- c(
    'Stability'
    , 'Validity'
    , 'Dist. (Classical)', 'Dist. (Opsahl)'
    , 'Discriminability'
)
colnames(test.mat) <- rep(test.names, times = 3)
# Set number of digits and replace undefined entries with empty strings
wh.na <- which(is.na(test.mat))
test.tab <- format(round(test.mat, 3), nsmall = 3)
test.tab[wh.na] <- ''
test.tab[1:length(test.tab)] <- paste0('\\(', test.tab, '\\)')

# Print table
file <- paste0(tabloc, 'tab-test.txt')
write(paste(
    '& ',
    paste(paste0('\\multicolumn{3}{c|}{', stat.incl, '}'), collapse = ' & '),
    ' \\\\'
), file = file)
write(paste(
    '& ',
    paste(colnames(test.mat), collapse = ' & '),
    ' \\\\\\hline'
), file = file, append = TRUE)
write.table(test.tab, file = file, append = TRUE,
            quote = FALSE, sep = ' & ', eol =  ' \\\\\n',
            row.names = TRUE, col.names = FALSE)

# Dynamic triadic closure of DDGGS1
dtc.ddggs1 <- dyn.triadic.closure(example$DDGGS1)
xtr.ddggs1 <- excl.transitivity(example$DDGGS1, type = 'global', stat = 'trans')
print(paste('Dynamic triadic closure of DDGGS1 =', round(dtc.ddggs1, 3)))
print(paste('Exclusive transitivity ratio of DDGGS1 =',
            round(xtr.ddggs1, 3)))

rm(list = ls())
