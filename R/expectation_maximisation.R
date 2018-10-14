#' EM steps

maximise <- function(Q,n,k,E){
  alpha <- sum(E*Q)/(k*sum(Q)) # division by 2 cancels out
  M = E*(1-Q)
  N = k*(1-Q)
  beta <- sum(M[lower.tri(M)])/ sum(N[lower.tri(N)])
  rho <- (1/choose(n,2))*sum(Q[lower.tri(Q)]) # false is the default
  rho <- (1/choose(n,2))*sum(Q[lower.tri(Q)]) # false is the default
  out <- c(alpha, beta, rho)
}

expectation <- function(alpha, beta, rho, E, k){
  Q = (rho*(alpha^E)*((1-alpha)^(k-E)))/(rho*(alpha^E)*((1-alpha)^(k-E))+(1-rho)*(beta^E)*((1-beta)^(k-E)))
}

EM <- function(alpha0, beta0, rho0, n, k, E){
  threshold=0.001
  alpha = alpha0
  beta = beta0
  rho = rho0
  prev=c(alpha, beta, rho)
  out=c(Inf,Inf,Inf)
  nIter=0
  while(!all(abs(prev-out)<c(threshold,threshold,threshold))){
    nIter=nIter+1
    prev <- c(alpha, beta, rho)
    Q <- expectation(alpha, beta, rho, E, k)
    out <- maximise(Q, n, k, E)
    alpha <- out[1]
    beta <- out[2]
    rho <- out[3]
  }
  out <- list(Q,alpha,beta,rho,nIter)
}





