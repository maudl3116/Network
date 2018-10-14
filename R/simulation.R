#' Simulating
#'
#' This script allows you to simulate network data and perform inference via the Expectation Maximization algorithm

n=50
k=5
p=0.1
output <- sampleErdosRenyi(n,p)
g <- output[[2]]
A <- output[[1]]
alpha <- 0.6 # true positive rate
beta <- 0.009 # false positive rate
E <- interact(A,alpha=0.6,beta=0.009, n,k)
E <- interact(A,alpha=0.6,beta=0.009, n,k)

simulation1 <- EM(alpha0=0.4, beta0=0.02, rho0=0.15, n, k, E)

