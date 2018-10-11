
install.packages("igraph")
library(igraph)

generate_graph <- function(n=50,p=0.1){
  g <- erdos.renyi.game(n, p, type = c("gnp"), directed=FALSE, loops=FALSE)
  plot(g, layout=layout.auto, vertex.size=0.8, vertex.label=NA, edge.arrow.size=0.2)
  A <- as.matrix(get.adjacency(g, type="lower"))
  A
}



interact <- function(A,alpha=0.6,beta=0.009, n=50,k){

  E = matrix(0,n,n)

  tmp = A[lower.tri(A)] # vector

  a = rbinom(n*(n-1)/2,k,alpha*tmp+beta*(1-tmp))

  E[lower.tri(E)] = a

  return(E)

}


maximise <- function(Q,n,k,E){
  alpha <- sum(E*Q)/(k*sum(Q)) # division by 2 cancels out
  M = E*(1-Q)
  N = k*(1-Q)
  beta <- sum(M[lower.tri(M)])/ sum(N[lower.tri(N)])
  rho <- (1/choose(n,2))*sum(Q[lower.tri(Q)]) # false is the default
  out <- c(alpha, beta, rho)
}

expectation <- function(alpha, beta, rho, E, k){
  Q = (rho*(alpha^E)*((1-alpha)^(k-E)))/(rho*(alpha^E)*((1-alpha)^(k-E))+(1-rho)*(beta^E)*((1-beta)^(k-E)))
}

EM <- function(alpha0, beta0, rho0, n, k, E){
  alpha = alpha0
  beta = beta0
  rho = rho0
  for (i in 1:100){
    Q <- expectation(alpha, beta, rho, E, k)
    out <- maximise(Q, n, k, E)
    alpha <- out[1]
    beta <- out[2]
    rho <- out[3]
  }
  print(alpha)
  print(beta)
  print(rho)
  out <- c(Q,alpha,beta,rho)
}

A <- generate_graph()
alpha <- 0.6 # true positive rate
beta <- 0.009 # false positive rate
E <- interact(A,alpha=0.6,beta=0.009, n=50,k=5)

simulation1 <- EM(alpha0=0.4, beta0=0.02, rho0=0.15, n=50, k=5, E)
simulation1

