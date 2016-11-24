# Global diagnostics of Mathematical Reviews network,
# within the aggregate, pure, and applied subnetworks,
# along a fixed-duration sliding window
library(Matrix)
library(igraph)
library(bitriad)
source('code/triadic-spec.R')
source('code/mathrev2igraph.R')
source('code/bonacich.centrality.R')

# Load MR data
load("calc/mathrev.RData")

pure2 <- sprintf('%02d', 3:59)
applied2 <- sprintf('%02d', 60:96)

# Aggregate MR centrality measures
mathrev.cent <- list()
for (yr in yrs) {
  bigraph <- as_an(paper.author.graph(mathrev[
    mathrev$year %in% (yr - dur + 1):yr, ]))
  # Degree centrality
  deg <- degree(actor_projection(bigraph))
  # Multidegree (2-walk) centrality
  mdeg.cent <- unit.cent(bigraph, 1, sparse = TRUE)
  # 4-walk centrality
  p.cent <- unit.cent(bigraph, 2, sparse = TRUE)
  # Eigenvector centrality
  eigen.cent <- ev.cent(bigraph, sparse = TRUE)
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
  mathrev.cent <- c(mathrev.cent, list(cent))
}

# Pure MR centrality measures
mathrev.pure <- mathrev[substr(mathrev$pclass, 1, 2) %in% pure2, ]
mathrev.pure.cent <- list()
for (yr in yrs) {
  bigraph <- as_an(paper.author.graph(mathrev.pure[
    mathrev.pure$year %in% (yr - dur + 1):yr, ]))
  # Degree centrality
  deg <- degree(actor_projection(bigraph))
  # Multidegree (2-walk) centrality
  mdeg.cent <- unit.cent(bigraph, 1, sparse = TRUE)
  # 4-walk centrality
  p.cent <- unit.cent(bigraph, 2, sparse = TRUE)
  # Eigenvector centrality
  eigen.cent <- ev.cent(bigraph, sparse = TRUE)
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
  mathrev.pure.cent <- c(mathrev.pure.cent, list(cent))
}

# Applied MR centrality measures
mathrev.applied <- mathrev[substr(mathrev$pclass, 1, 2) %in% applied2, ]
mathrev.applied.cent <- list()
for (yr in yrs) {
  bigraph <- as_an(paper.author.graph(mathrev.applied[
    mathrev.applied$year %in% (yr - dur + 1):yr, ]))
  # Degree centrality
  deg <- degree(actor_projection(bigraph))
  # Multidegree (2-walk) centrality
  mdeg.cent <- unit.cent(bigraph, 1, sparse = TRUE)
  # 4-walk centrality
  p.cent <- unit.cent(bigraph, 2, sparse = TRUE)
  # Eigenvector centrality
  eigen.cent <- ev.cent(bigraph, sparse = TRUE)
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
  mathrev.applied.cent <- c(mathrev.applied.cent, list(cent))
}

save(mathrev.cent, mathrev.pure.cent, mathrev.applied.cent,
     file = 'calc/mathrev-center.RData')

rm(list = ls())
