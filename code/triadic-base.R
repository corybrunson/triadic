# Preliminaries for tmtc image generation
source('code/triadic-spec.R')

# Destination
figloc <- 'fig/'
tabloc <- 'tab/'

# Settings
eps <- FALSE
options(scipen = 10)
wid <- 4

# Image-generating shell
img <- if(eps) postscript else pdf
suf <- if(eps) '.eps' else '.pdf'

# Graph aesthetics
g.shape = c('circle', 'square')
g.color = c('SkyBlue2', 'lightcoral')
g.bw = c('gray85', 'gray75')
g.size = c(28, 24)
edge.wid <- 2

# Triangular numbers
triangular <- function(n) n * (n + 1) / 2
undo.triangular <- function(N) (-1 + sqrt(1 + 8 * N)) / 2

# Fonts for axes
quartzFonts(sans = quartzFont(c('HelveticaNeue', 'HelveticaNeue-Bold',
                                'HelveticaNeue-Italic', 'HelveticaNeue-Bold')))

# Wobble layout coordinates a bit
library(MASS)
wobble.coords <- function(l, S) l + t(sapply(1:dim(l)[1], function(x) {
    mvrnorm(n = dim(l)[2], mu = 0, Sigma = S)
}))

# Convertinge wedge lists to clustering coefficients
wedges2cc <- function(wl, v, t) sum(wl[, t]) / sum(wl[, v])
# Converting wedge lists to transitivity ratios
wedges2tr <- function(wl, v, t) sum(wl[, t] / (3 * sum(wl[, v])))

# Asterisks according to p-value
p2ast <- function(p) if(p >= .01) '' else
    paste('{\\small',
          paste(rep('*', min(4, ceiling(-log(p, 10)) - 2)), collapse = ''),
          '}',
          sep = '')

# TABLE SHELLS

# LaTeX tables for *Network Science* specifications
latex.table <- function(x, digits, align, file, math.mode = c(),
                        hline.extra = c(), ...) {
    x <- as.matrix(x)
    if(!is.null(math.mode)) {
        if(any(grepl('all|entries', math.mode)))
            x[1:length(x)] <- paste0('\\(', x, '\\)')
        if(any(grepl('all|row|^names', math.mode)))
            rownames(x) <- paste0('\\(', rownames(x), '\\)')
        if(any(grepl('all|col|^names', math.mode)))
            colnames(x) <- paste0('\\(', colnames(x), '\\)')
    }
    print(xtable(x = x, digits = digits, align = align), file = file,
          floating = FALSE,
          hline.after = c(-1, -1, 0, hline.extra, nrow(x), nrow(x)),
          sanitize.text.function = identity, ...)
}

latex.tables <- function(
    list, digits, align.row, align.each, align.headers, headerline = FALSE,
    file, math.mode = c(), hline.extra = c(), ...) {
    x <- do.call(cbind, list)
    if(!is.null(math.mode)) {
        if(any(grepl('all|entries', math.mode)))
            x[1:length(x)] <- paste0('\\(', x, '\\)')
        if(any(grepl('all|row|^names', math.mode)))
            rownames(x) <- paste0('\\(', rownames(x), '\\)')
        if(any(grepl('all|col|^names', math.mode)))
            colnames(x) <- paste0('\\(', colnames(x), '\\)')
    }
    list.line <- paste(paste(c('\\hline\\hline \n',
                               paste0('\\multicolumn{',
                                      sapply(list, ncol), '}{',
                                      align.headers, '}{',
                                      names(list), '}')), collapse = ' & '),
                       '\\\\', ifelse(headerline, '\\hline', ''), ' \n')
    # Print table to file with headers line appended to a row
    print(xtable(x = x, digits = digits,
                 align = paste0(c(align.row, rep(align.each, length(list))),
                                collapse = '')), file = file,
          hline.after = c(0, hline.extra, nrow(x), nrow(x)),
          floating = FALSE, sanitize.text.function = identity, ...,
          add.to.row = list(list(-1), list.line))
}

# Use expressions for italic(C) and italic(TC) but strings elsewhere
label_custom <- function(variable, value) {
    ifelse(grepl('italic\\([CT]\\)', value),
           label_parsed(variable, value),
           label_value(variable, value))
}

# Turn a null model list of data frames into a grid of histograms
# (lst items must be named, and each must be a data frame with the same cols)
hist.grid <- function(lst, col = 'red', lty = 'dashed') {
    
    # Bind null data frames into a single data frame
    null.df <- do.call(
        rbind,
        lapply(1:length(lst), function(j) {
            cbind(indiv = names(lst)[j],
                  lst[[j]][-1, ],
                  type = 'null')
        }))
    melt.null <- melt(null.df, id.vars = c('indiv', 'W', 'type'))
    
    # Bind observations into a single data frame
    obs.df <- do.call(
        rbind,
        lapply(1:length(lst), function(j) {
            cbind(indiv = names(lst)[j],
                  lst[[j]][1, ],
                  type = 'obs')
        }))
    melt.obs <- melt(obs.df, id.vars = c('indiv', 'W', 'type'))
    
    # Construct plot
    plot <- ggplot(melt.null, aes(x = value)) +
        geom_histogram(aes(weight = W), position = 'identity',
                       binwidth = .1 / (log2(1000) + 1)) +
        geom_vline(data = melt.obs, aes(xintercept = value),
                   color = col, linetype = lty, size = 1) +
        facet_grid(indiv ~ variable, labeller = label_custom) +
        xlab('') + ylab('') +
        scale_x_continuous(labels = c(' 0', .25, .5, .75, ''),
                           expand = c(0, 0)) +
        scale_y_continuous(breaks = NULL, expand = c(0, 0)) +
        coord_cartesian(xlim = c(0, 1)) +
        theme_bw()
    print(plot)
    
}
