#' This function makes it possible to fit an ergm model (i.e. sample from the posterior distribution on the parameters) using an enhanced MCMC technique.
#' @param yobs one observation of the network
#' @param initial_theta The initial guess on the values of the parameters of the ergm model
#' @param lower_theta Used to bound the possible values of the parameters
#' @param upper_theta Used to bound the possible values of the parameters
#' @param mcmc_iter The number of iterations of the MCMC algorithm
#' @examples

exchange <- function(X, COV, mcmc_iter, thetas, a){

   n = ncol(X)

  statistics_x <- Summary(X)

  old_theta <- thetas[1,]


  for(iter in seq(2:mcmc_iter+1)){

    # sample new parameters to propose
    proposed_theta <- rmvnorm(n=1, mean=old_theta, COV, method="chol")

    proposed_theta <- as.vector(proposed_theta)

    # sample auxiliary variable from the likelihood at theta=proposed_theta
    proposed_u <- simulate(network(n,directed=FALSE)~edges+kstar(2) + kstar(3)+triangles, coef=proposed_theta)
    proposed_u <- as.matrix(proposed_u)
    # compute h()
    statistics_proposed_u <- Summary(proposed_u)
    log_ratio_u <- sum(old_theta*statistics_proposed_u)-sum(proposed_theta*statistics_proposed_u)
    log_ratio_x <- sum(proposed_theta*statistics_x)-sum(old_theta*statistics_x)

    # compute the other values appearing in the log acceptance ratio

    # priors
    log_prior_theta_old <- log_prior_theta(old_theta,a)
    log_prior_theta_proposed <- log_prior_theta(proposed_theta,a)


    # form the acceptance rate
    r <-log_prior_theta_proposed - log_prior_theta_old + log_ratio_u + log_ratio_x

      # accept-reject
    tmp = runif(1)

    if( log(tmp) < r){
        old_theta <- proposed_theta
    }
    else{
      old_theta <- old_theta
    }


    thetas[iter,] = old_theta

    if (iter%%100==0){
        print(iter)
        remove = c(iter:length(thetas[,1]))
        plot(thetas[-remove,1],type="l")
        plot(thetas[-remove,2],type="l")
      }
  } # end iterations
  return(thetas)
}



log_prior_theta <- function(theta,a){
  theta_1 <- dunif(theta[1],-a,a,log=TRUE)
  theta_2 <- dunif(theta[2],-a,a,log=TRUE)
  theta_3 <- dunif(theta[3],-a,a,log=TRUE)
  theta_4 <- dunif(theta[4],-a,a,log=TRUE)
  return(theta_1+theta_2+theta_3+theta_4)
}

Summary = function(A){
  n = nrow(A)
  star = apply(A,1,sum)
  result = numeric(4)

  for(i in seq(1:n)){
    result[1] = result[1] + choose(star[i],1)
    result[2] = result[2] + choose(star[i],2)
    result[3] = result[3] + choose(star[i],3)
  }
  result[1]=result[1]/2
  for(k in seq(1:n)){
    for(j in seq(1:k+1)){
      for(i in seq(1:j+1)){
        result[4] = result[4] + A[i,k]*A[k,j]*A[j,i] }}}

  return(result)
}

