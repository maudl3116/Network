
convergence <- function (N=200, alpha0=0.5, beta0=0.05, rho0=0.2, n=50, k=10, E){
  E <- simulate_infer()[[3]]
  matrix=matrix(0,N,3)
  alpha = alpha0
  beta = beta0
  rho = rho0
  for (i in 1:N){
    Q <- expectation(alpha, beta, rho, E, k)
    out <- maximise(Q, n, k, E)
    alpha <- out[1]
    beta <- out[2]
    rho <- out[3]
    matrix[i,]=c(alpha, beta, rho)
    }
  m1<-matrix[-c(1:50),]
  parmfrow=c(3,1)
  plot(m1[,1], type="l")
  plot(m1[,2], type="l")
  plot(m1[,3], type="l")
  parmfrow=c(1,1)
  return(m1[N,])
}
