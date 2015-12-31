# Bonacich centrality scores

# Compute unit vector for centrality due to paths of length 2 * k
unit.cent <- function(bigraph, k, beta, rm.diag = FALSE, sparse = FALSE) {
    stopifnot(k >= 1)
    # Incidence matrix
    inc <- as_incidence_matrix(bigraph, sparse = sparse)
    # Weighted adjacency matrix (transpose cross product)
    mat <- if (!sparse) {
        inc %*% t(inc)
    } else {
        tcrossprod(inc)
    }
    if(rm.diag) mat <- mat - diag(diag(mat))
    # Default: beta approximately equal to first eigenvalue of mat
    if(missing(beta)) {
        lambda <- if (!sparse) {
            eigen(mat)$values[1]
        } else {
            lambda <- arpack(function(x, extra = NULL) as.vector(mat %*% x),
                             options = list(n = nrow(mat),
                                            nev = 2,
                                            ncv = 5,
                                            which = "LM",
                                            maxiter = 500),
                             sym = TRUE, complex = FALSE)$values[1]
        }
        if(any(lambda < 0)) lambda <- -lambda
        beta <- 100000 / ceiling(100000 * lambda)
    }
    # Construct mat ^ i * beta ^ (i - 1) for i = 1, ..., k
    mats <- list()
    i <- 1; pow <- mat
    mats <- c(mats, list(pow))
    while(i < k) {
        i <- i + 1
        pow <- beta * (pow %*% mat)
        mats <- c(mats, list(pow))
    }
    # Add these powers of mat
    bigmat <- Reduce('+', mats)
    # Take the row sums
    cent <- rowSums(bigmat)
    # Normalize so that the length of cent is 1
    cent <- cent / min(cent) # necessary to avoid rounding errors
    cent / sqrt(sum(cent ^ 2))
}

# Compute eigenvector centrality
ev.cent <- function(bigraph, rm.diag = FALSE, sparse = FALSE) {
    # Incidence matrix
    inc <- as_incidence_matrix(bigraph, sparse = sparse)
    # Weighted adjacency matrix
    mat <- inc %*% t(inc)
    if(rm.diag) mat <- mat - diag(diag(mat))
    # Vector of eigenvector centralities
    eigen.cent <- if (!sparse) {
        eigen(mat)$vector[, 1]
    } else {
        arpack(function(x, extra = NULL) as.vector(mat %*% x),
               options = list(n = nrow(mat),
                              nev = 1,
                              ncv = 5,
                              which = "LM",
                              maxiter = 500),
               sym = TRUE, complex = FALSE)$vector
    }
    if(any(eigen.cent < -10^-10)) eigen.cent <- -eigen.cent
    stopifnot(all(eigen.cent >= -10^-10))
    stopifnot(abs(sum(eigen.cent ^ 2) - 1) < 10^-10)
    eigen.cent
}
