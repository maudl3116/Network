#' A Data Generator
#'
#' This script contains different functions to generate graphs from different Exponential Random Graph Models


if(!require(igraph)){
  install.packages("somepackage")
  library(igraph)
}

#' This function generates an Erdos Renyi graph, the model used in the paper
#' @param n The number of nodes in the network
#' @param p The probability of having an edge at any position (i,j)
#' @examples
#' sampleErdosRenyi_function()
sampleErdosRenyi <- function(n=50,p=0.1){
  g <- erdos.renyi.game(n, p, type = c("gnp"), directed=FALSE, loops=FALSE)
  plot(g, layout=layout.auto, vertex.size=6, vertex.label=NA, edge.arrow.size=0.2,vertex.color="light blue")
  A <- as.matrix(get.adjacency(g, type="lower"))
  list(A,g)
}

#' This function generates a graph from an ERGM
#' @param n The number of nodes in the network

ERGM <- function(n, param){
  # network to initialise the MCMC algorithm
  g0 <- network(n,directed=FALSE)
  g <- simulate(~edges+triangles, nsim=1, coef=param, basis=g0, control=control.simulate(MCMC.burnin=100,MCMC.interval=100))
  plot(g, vertex.size=6, vertex.label=NA, edge.arrow.size=0.2,vertex.color="light blue")
}

#' This function generates noisy data
#' @param A The underlying ERGM
#' @param alpha The true positive rate, i.e. the probability of observing an edge at (i,j) when A_{(i,j)}=1
#' @param beta The false positive rate, i.e. the probability of observing an edge at (i,j) when A_{(i,j)}=0
#' @param n The number of nodes in the network
#' @param k The number of observations for each pair (i,j)
#' @examples
interact <- function(A,alpha=0.6,beta=0.009, n=50,k){

  E = matrix(0,n,n)

  tmp = A[lower.tri(A)] # vector

  a = rbinom(n*(n-1)/2,k,alpha*tmp+beta*(1-tmp))

  E[lower.tri(E)] = a

  return(E)

}
