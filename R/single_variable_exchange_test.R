#' evaluating the exact-approximate MCMC algorithm on a concrete example for which all computations are easy
#' we want to estimate the precision parameter of a zero-mean univariate gaussian, based on one observation
#'
#' @param yobs one observation from a normal distribution
#' @param initial_theta The initial guess on the value of the precision of the normal
#' @param mcmc_iter The number of iterations of the exact-approximate MCMC algorithm
#' @param step The standard deviation used in the random walk


benchmark <- function(yobs, initial_theta, mcmc_iter, step=0.5){
  accepted=0
  old_theta <- initial_theta
  path = numeric(mcmc_iter)

  for(iter in seq(1:mcmc_iter)){
    # sample new parameters to propose
     proposed_theta <- rnorm(n=1, mean=old_theta, sd=step)
    # proposed_theta <- rgamma(1,0.5+1,0.5*yobs^2+1)  # sample from the posterior ! to benchmark

    if(proposed_theta <= 0){

      path[iter]=old_theta

    }else{


      # sample auxiliary variable from the likelihood at theta=proposed_theta
      proposed_u <- rnorm(1,mean=0,sd=sqrt(1/proposed_theta))
      # compute the approximation of the log-ratio of the partition functions
      log_ratio_partitions_est <- 0.5*(proposed_u^2)*(proposed_theta-old_theta)


      # compute the other values appearing in the log acceptance ratio

      # priors
      log_prior_theta_old <- dgamma(old_theta,1,1,log=TRUE)
      log_prior_theta_proposed <- dgamma(proposed_theta,1,1,log=TRUE)

      # transitions  (assume symmatric for now)
      #log_theta_old <- log_theta_kernel(proposed_theta, old_theta)
      #log_theta_proposed <- log_theta_kernel(old_theta, proposed_theta)
      #log_theta_old <- dgamma(0.5+1,0.5*yobs^2+1,log=TRUE)
      #log_theta_proposed <- dgamma(0.5+1,0.5*yobs^2+1,log=TRUE)
      log_theta_old <- 0
      log_theta_proposed <- 0

      # unnormalised likelihoods
      log_likelihood_ratio = 0.5*(yobs^2)*(old_theta-proposed_theta)


      # form the acceptance rate
      a <- log_likelihood_ratio + log_ratio_partitions_est + log_prior_theta_proposed - log_prior_theta_old + log_theta_proposed - log_theta_old

      # accept-reject
      tmp = runif(1)

      if( log(tmp) < a){
        old_theta <- proposed_theta
        if(iter>5000) accepted=accepted+1
      }


      path[iter]=proposed_theta
      if (iter%%1000==0){
        print(iter)
        remove = c(iter:length(path))
        plot(path[-remove],type="l")
      }
    }

  } # end iterations
  print(accepted)
  return(list(parameters=old_theta,estimates=path))
}

# run

gt_theta=0.5
yobs = 1

out=benchmark(yobs, initial_theta=1, mcmc_iter=10000,step=0.1)

print(mean(out$estimates[seq(9000,10000)]))
