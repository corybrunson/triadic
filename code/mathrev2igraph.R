# Convert Mathematical Reviews data to formats readable into igraph

# List of authors by paper
authors <- function(data) strsplit(data[, 3], '\\|')

# List of subject classifications by paper (NECESSARY?)
subjects <- function(data) strsplit(paste(
  data[, 4],
  data[, 5],
  sep = '|'
), '\\|')

# Bipartite paper-author graph (when directed, papers point to authors)
# (Paper nodes are labeled by primary classification and year of publication)
paper.author.graph <- function(data, directed = FALSE) {
  if (dim(data)[1] == 0) return(graph.empty(n = 0))
  data.authors <- authors(data)
  uniq.authors <- unique(unlist(data.authors))
  edf <- as.data.frame(matrix(unlist(sapply(mapply(
    paste,
    data$index,
    data.authors,
    sep = ','
  ), strsplit, split = ',')), nc = 2, byrow = TRUE))
  vdf <- data.frame(id = c(data$index, uniq.authors),
                    pclass = c(data$pclass, rep(NA, length(uniq.authors))),
                    time = c(data$year, rep(NA, length(uniq.authors))))
  graph <- graph.data.frame(edf, directed = directed, vertices = vdf)
  V(graph)$type <- !(substr(V(graph)$name, 3, 3) %in% letters)
  graph
}

# Bipartite paper-subject graph (when directed, papers point to subjects)
paper.subject.graph <- function(data, chars = 5, directed = FALSE) {
  if (dim(data)[1] == 0) return(graph.empty(n = 0))
  data.subjects <- lapply(subjects(data), substr, 1, chars)
  uniq.subjects <- unique(unlist(data.subjects))
  edf <- as.data.frame(matrix(unlist(sapply(mapply(
    paste,
    data$index,
    data.subjects,
    sep = ','
  ), strsplit, split = ',')), nc = 2, byrow = TRUE))
  vdf <- data.frame(id = c(data$index, uniq.subjects))
  graph <- graph.data.frame(edf, directed = directed, vertices = vdf)
  V(graph)$type <- nchar(V(graph)$name) > 5
  graph
}

# Unipartite author graph
author.graph <- function(
  data, weighted = FALSE, wt = function(vec) 1, wtfun = sum, multiedge = FALSE
) {
  if(dim(data)[1] == 0) return(graph.empty(n = 0))
  alst <- authors(data)
  if(max(sapply(alst, length)) == 1) {
    uniq <- unique(unlist(alst))
    graph <- graph.empty(n = length(uniq))
    V(graph)$name <- uniq
    return(graph)
  }
  el <- as.data.frame(t(matrix(unlist(lapply(
    alst,
    function(vec) if(length(vec) == 1) c() else combn(vec, m = 2)
  )), nr = 2)))
  if(multiedge) el <- cbind(el, year = rep(data$year, sapply(
    alst,
    function(vec) if(length(vec) == 1) 0 else choose(length(vec), 2)
  )))
  graph <- graph.data.frame(el, directed = FALSE,
                            vertices = data.frame(unique(unlist(alst))))
  if(weighted) {
    E(graph)$weight <- rep(sapply(alst, wt),
                           choose(sapply(alst, length), 2))
    simplify(graph, edge.attr.comb = list(weight = wtfun))
  } else if(multiedge) graph else simplify(graph)
}

# Unipartite subject graph
subject.graph <- function(
  data, chars = 5, weighted = FALSE, wt = function(vec) 1, wtfun = sum
) {
  if (dim(data)[1] == 0) return(graph.empty(n = 0))
  slst <- lapply(subjects(data), substr, 1, chars)
  el <- t(matrix(unlist(lapply(
    slst,
    function(vec) if(length(vec) == 1) c() else combn(vec, m = 2)
  )), nr = 2))
  graph <- graph.edgelist(el, directed = FALSE)
  if(weighted) {
    E(graph)$weight <- rep(sapply(alst, wt),
                           choose(sapply(alst, length), 2))
    simplify(graph, edge.attr.comb = list(weight = wtfun))
  } else simplify(graph)
}

# Data set contributing to largest component
lc.data <- function(data) {
  data[sapply(
    data$authors,
    substr, start = 1, stop = 9
  ) %in% V(largest.component(author.graph(data)))$name, ]
}

# Sort character strings into pure and applied
pure.or.applied <- function(string) {
  sapply(as.numeric(substr(string, 1, 2)), function(x) if(is.na(x)) NA else
    if(x > 2 & x < 60) 'pure' else if(x > 59 & x < 97) 'applied' else NA)
}

# Colored unipartite publication graph
paper.graph <- function(data, char = 'pureapplied',
                        multiplicity = TRUE, keep = TRUE) {
  bigraph <- paper.author.graph(data)
  stopifnot(all(data$index == V(bigraph)$name[V(bigraph)$type]))
  V(bigraph)$sub[V(bigraph)$type] <- if(char == 'pureapplied')
    pure.or.applied(data$pclass) else as.numeric(substr(data$pclass, 1, char))
  graph <- bipartite.projection(bigraph, multiplicity = multiplicity)[[1]]
  if(!keep) graph <- delete.vertices(graph, which(is.na(V(graph)$sub)))
  graph
}

# productivity assortativity
assortativity.productivity <- function(data, baseline = FALSE, rm1 = FALSE) {
  ppa <- table(unlist(authors(data)))
  wgraph <- author.graph(data, weighted = TRUE)
  ppa <- ppa[unlist(sapply(1:length(ppa), function(i) {
    which(names(ppa) == V(wgraph)$name[i])
  }))]
  if(baseline) return(assortativity(wgraph, types1 = ppa))
  m <- ecount(wgraph)
  edges <- get.edgelist(wgraph, names = FALSE)
  cwt <- E(wgraph)$weight
  if(rm1) {
    ww <- which((ppa[edges[, 1]] > 1) & (ppa[edges[, 2]] > 1))
    edges <- edges[ww, ]
    cwt <- cwt[ww]
  }
  rem1 <- ppa[edges[, 1]] - cwt
  rem2 <- ppa[edges[, 2]] - cwt
  num1 <- sum(rem1 * rem2) / m
  num2 <- (sum(rem1 + rem2) / (2 * m)) ^ 2
  dem <- sum(rem1 ^ 2 + rem2 ^ 2) / (2 * m)
  return((num1 + num2) / (dem + num2))
}

# FUNCTION: Given MR data for window + increment, compute two vectors:
# L = vector indexed by node of number of pairs of unlinked neighbors in G_w
# C = vector indexed by node of number of these pairs linked in G_(w+i)
# Global proportion is sum(C) / sum(L)
# Average opportunity-dependent local proportion is
# sum(C[degree(g) == k]) / sum(L[degree(g) == k])
author.edgelist <- function(data) {
    au <- unlist(lapply(authors(data), function(vec) {
        if(length(vec) == 1) c() else combn(vec, m = 2)
    }))
    if(is.null(unlist(unique(au)))) return(matrix('', nr = 0, nc = 2))
    return(as.data.frame(t(matrix(au, nr = 2)), stringsAsFactors = FALSE))
}

# PROBLEM: denom is number of pairs of distance 2 in g.0, not number of wedges
vee.closure <- function(data, window, increment) {
    el.0 <- author.edgelist(data[data$year %in% window, ])
    if(dim(el.0)[1] == 0) return(NA)
    el.1 <- author.edgelist(data[data$year %in% increment, ])
    g.0 <- simplify(graph.data.frame(el.0, directed = FALSE))
    denom <- sum(sapply(neighborhood(g.0, order = 2), length) -
                     degree(g.0) - 1) / 2
    el <- data.frame(V1 = c(el.0[, 1], el.1[, 1]),
                     V2 = c(el.0[, 2], el.1[, 2]),
                     old = c(rep(TRUE, dim(el.0)[1]), rep(FALSE, dim(el.1)[1])))
    g <- simplify(graph.data.frame(el, directed = FALSE),
                  edge.attr.comb = list(old = any))
    E.new <- which(!E(g)$old)
    num <- sum(sapply(E.new, function(e) {
        ge <- get.edge(g, e)
        if(!(all(ge <= vcount(g.0)))) 0 else
            (shortest.paths(g.0, ge[1], ge[2]) == 2)
    }))
    num / denom
}

# Alternative: check neighborhoods rather than distance, include local option
wedge.closure <- function(data, window, increment, type = 'global') {
    # Two-column list of edges in "present" graph using author names
    el.0 <- author.edgelist(data[data$year %in% window, ])
    if(dim(el.0)[1] == 0) return(NA)
    # Two-column list of edges in "increment" graph using author names
    el.1 <- author.edgelist(data[data$year %in% increment, ])
    # Graph object based on "present" edges
    g.0 <- simplify(graph.data.frame(el.0, directed = FALSE))
    # The denominator is the number of open wedges at each node
    denom <- round(transitivity(g.0, type = 'local') * choose(degree(g.0), 2))
    
    # Combined edge list with indicator for edges in g.0
    el <- data.frame(V1 = c(el.0[, 1], el.1[, 1]),
                     V2 = c(el.0[, 2], el.1[, 2]),
                     old = c(rep(TRUE, dim(el.0)[1]), rep(FALSE, dim(el.1)[1])))
    # Graph object based on all edges with indicator for edges extant in g.0
    g <- simplify(graph.data.frame(el, directed = FALSE),
                  edge.attr.comb = list(old = any))
    # Vector of nodes of g.0 (if any) corresponding to nodes of g
    v.0 <- match(V(g)$name, V(g.0)$name)
    # Identify nodes with wedges in g.0 closed by each new edge in g
    vs.0 <- lapply(which(!E(g)$old), function(e) {
        ge <- get.edges(g, e)
        if(any(is.na(v.0[ge]))) return(c()) else
            do.call(intersect, neighborhood(g.0, v.0[ge], order = 1))
    })
    # The nuerator will be the number of newly closed wedges at each node of g.0
    num <- if(is.null(unlist(vs.0))) rep(0, vcount(g.0)) else
        tabulate(unlist(vs.0), nbins = vcount(g.0))
    
    wh <- which(is.na(denom))
    stopifnot(sum(num[wh]) == 0)
    denom[wh] <- 0
    
    # Mimic twomode.transitivity output
    wedges <- rbind(denom, num)
    if(type == 'global') return(sum(wedges[2, ]) / sum(wedges[1, ]))
    if(type == 'local') return(wedges[2, ] / wedges[1, ])
    return(data.frame(V = wedges[1, ], T = wedges[2, ]))
}
