#' This function makes it possible to fit an ergm model (i.e. sample from the posterior distribution on the parameters) using an enhanced MCMC technique.
#' @param yobs one observation of the network
#' @param initial_theta The initial guess on the values of the parameters of the ergm model
#' @param lower_theta Used to bound the possible values of the parameters
#' @param upper_theta Used to bound the possible values of the parameters
#' @param mcmc_iter The number of iterations of the MCMC algorithm
#' @examples

ergm_no_noise_fit <- function(n, yobs, initial_theta, lower_theta, upper_theta, mcmc_iter){

  # to do:
  # test with different priors
  # put several observations


  old_theta <- initial_theta

  m <- matrix(rbinom(n*n,1,.5),n,n)
  diag(m) <- 0
  net_0=network(m,directed=FALSE)
  path = matrix(0,mcmc_iter,2)

  for(iter in seq(1:mcmc_iter)){

    # sample new parameters to propose
    proposed_theta <- trunc_norm(old_theta, lower_theta, upper_theta)
    proposed_theta <- as.vector(proposed_theta)
    # sample auxiliary variable from the likelihood at theta=proposed_theta
    proposed_u <- simulate(~edges+triangles,coef=proposed_theta,basis=net_0,control=control.simulate(MCMC.burnin=100,MCMC.interval=100))

    # compute the approximation of the log-ratio of the partition functions
    statistics_proposed_u <- statistics_network(proposed_u)
    log_ratio_partitions_est <- sum(old_theta*statistics_proposed_u)-sum(proposed_theta*statistics_proposed_u)

    # compute the other values appearing in the log acceptance ratio

      # priors
    log_prior_theta_old <- log_prior_theta(old_theta)
    log_prior_theta_proposed <- log_prior_theta(proposed_theta)

      # transitions  (assume symmatric for now)
    log_theta_old <- log_theta_kernel(proposed_theta, old_theta)
    log_theta_proposed <- log_theta_kernel(old_theta, proposed_theta)
    #log_theta_old <- 0
    #log_theta_proposed <- 0


      # unnormalised likelihoods
    yobs_stats <- statistics_network(yobs)
    log_likelihood_old <- sum(yobs_stats*old_theta)
    log_likelihood_proposed <- sum(yobs_stats*proposed_theta)


      # form the acceptance rate
    v <- c(0,log_likelihood_proposed - log_likelihood_old + log_ratio_partitions_est + log_prior_theta_proposed - log_prior_theta_old + log_theta_proposed - log_theta_old)

      # accept-reject
    tmp = runif(1)

    if( log(tmp) < min(v) ){
      old_theta <- proposed_theta
    }


    path[iter,]=old_theta
    if (iter%%100==0){
      print(iter)
      remove = c(iter:length(path[,1]))
      plot(path[-remove,1],type="l")
      plot(path[-remove,2],type="l")
    }

  } # end iterations
  return(list(parameters=old_theta,estimates=path))
}

trunc_norm <- function(old_theta, lower_theta, upper_theta){
  #out <- rtmvnorm(n=1, mean=old_theta, lower=lower_theta, upper=upper_theta, algorithm="rejection")
  out <-  rmvnorm(n=1, mean=old_theta, sigma=2*diag(length(old_theta)), method="chol")
  }

statistics_network <- function(proposed_u){

 g <- asIgraph(proposed_u)
 S1 <- network.edgecount(proposed_u)
 S2 <- length(cliques(g,min=3,max=3))

}

log_prior_theta <- function(theta){
 #theta_1 <- dnorm(theta[1], mean=0, sd=sqrt(30), log=TRUE)
 #theta_2 <- dnorm(theta[2], mean=0, sd=sqrt(30), log=TRUE)
theta_1 <- dunif(theta[1],-8,0,log=TRUE)
theta_2 <- dunif(theta[2],0,10,log=TRUE)
 return(theta_1+theta_2)
}

log_theta_kernel <- function(proposed_theta, old_theta){
  dmvnorm(proposed_theta, mean=old_theta, sigma=2*diag(length(old_theta)), log=TRUE)
}

