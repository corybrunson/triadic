# Expected clustering coefficient on bipartite random graph model
moment <- function(dist, n, mode = 'frequency') {
    vec <- as.numeric(switch(mode,
                             'frequency' = dist,
                             'value' = tabulate(dist)))
    sum(((1:length(vec)) ^ n) * vec)
}
nsw.clust <- function(actor.degree, event.degree) {
    1 / (((moment(actor.degree, 2) - moment(actor.degree, 1)) *
              (moment(event.degree, 2) - moment(event.degree, 1)) ^ 2) /
             (moment(actor.degree, 1) * moment(event.degree, 1) *
                  (2 * moment(event.degree, 1) - 3 * moment(event.degree, 2) +
                       moment(event.degree, 3))) + 1)
}

# Use only first actor in each structural equivalence class
local.reps <- function(bigraph) {
    mat <- get.incidence(bigraph)
    reps <- setdiff(1:nrow(mat),
                    which(sapply(2:nrow(mat), function(j) {
                        any(sapply(1:(j - 1), function(k) {
                            all(mat[j, ] == mat[k, ])
                        }))})) + 1)
    names(reps) <- rownames(mat)[reps]
    return(reps)
}

# Convert matrices to network objects, perform SISs, export results as matrices
mat.sis <- function(mat, nsim = 1, save.mats = TRUE, save.labels = FALSE) {
    
    # Use tall, thin matrices
    mat.fn <- if(ncol(mat) > nrow(mat)) t else identity
    
    # Build network object (from tall and thin matrix)
    net <- network(mat.fn(mat))
    net %v% 'set' <- c(rep(1, nrow(mat)), rep(2, ncol(mat)))
    
    # Run simulations and store network objects
    sim <- simulate(net, nsim = nsim, save.networks = save.mats)
    
    # Convert network objects to matrices
    sim[[1]] <- if(is.network(sim[[1]])) mat.fn(as.matrix(sim[[1]])) else
        lapply(sim[[1]], function(s) mat.fn(as.matrix(s)))
    if(save.labels) if(is.matrix(sim[[1]])) {
        rownames(sim[[1]]) <- rownames(mat)
        colnames(sim[[1]]) <- colnames(mat)
    } else
        for(i in 1:length(sim[[1]])) {
            rownames(sim[[1]][[i]]) <- rownames(mat)
            colnames(sim[[1]][[i]]) <- colnames(mat)
        }
    return(sim)
}

# Perform SIS on igraph objects; save all produced
igraph.sis <- function(bigraph, nsim = 1, save.labels = FALSE) {
    mat <- get.incidence(bigraph)
    sim <- mat.sis(mat, nsim = nsim,
                   save.mats = TRUE, save.labels = save.labels)
    sim[[1]] <- if(is.matrix(sim[[1]])) graph.incidence(sim[[1]]) else
        lapply(sim[[1]], graph.incidence)
    return(sim)
}

if(FALSE) {
    
    # 1. Density model
    
    random.bigraph.density <- function(n1, n2, m, type = 'm', size = 1) {
        library(igraph)
        if(type == 'p') {
            mats <- lapply(1:size, function(i)
                matrix(replicate(n1 * n2, rbinom(1, 1, m / (n1 * n2)))))
        } else {
            mats <- lapply(1:size, function(i) matrix(0, nr = n1, nc = n2))
            for(i in 1:size) mats[[i]][sample(1:(n1 * n2), m)] <- 1
        }
        return(lapply(mats, graph.incidence))
    }
    
    # 2. Standard configuration model
    
    mat.sis.sample <- function(rowsum, colsum, nsim = 1, save.labels = FALSE) {
        library(networksis)
        
        # Build network object
        net <- network(mat)
        net %v% 'set' <- c(rep(1, nrow(mat)), rep(2, ncol(mat)))
        
        # Run simulations and store network objects
        sim <- simulate(net, nsim = nsim, save.networks = save.mats)[[1]]
        
        # Convert network objects to matrices
        sim <- if(is.network(sim)) {
            list(as.matrix(sim))
        } else lapply(sim, as.matrix)
        if(save.labels) if(is.matrix(sim)) {
            rownames(sim) <- rownames(mat)
            colnames(sim) <- colnames(mat)
        } else
            for(i in 1:length(sim)) {
                rownames(sim[[i]]) <- rownames(mat)
                colnames(sim[[i]]) <- colnames(mat)
            }
        return(sim)
    }
    
    random.bigraph.config <- function(deg1, deg2, type = 'm', size = 1) {
        if(type == 'p') {
            mats <- mat.sis.sample(deg1, deg2, nsim = size)
        } else {
            stop('wiring/subset option not yet implemented')
        }
        return(lapply(mats, graph.incidence))
    }
    
    # 3. Degree-degree profile-preserving configuration model
    
    mat.sis.profile.sample <- function(mat, nsim = 1, save.labels = FALSE) {
        
        # Across pairs of degrees with more than one representative each...
        
        # Conduct a SIS of their submatrix
        
    }
    
    random.bigraph.profile <- function(mat, type = 'm', size = 1) {
        if(type == 'p') {
            mats <- mat.sis.profile.sample(mat, nsim = size)
        } else {
            stop('wiring/subset option not yet implemented')
        }
        return(lapply(mats, graph.incidence))
    }
    
}
