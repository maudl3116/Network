#' This function makes it possible to fit an ergm model (i.e. sample from the posterior distribution on the parameters) using an enhanced MCMC technique.
#' @param yobs one observation of the network
#' @param initial_theta The initial guess on the values of the parameters of the ergm model
#' @param lower_theta Used to bound the possible values of the parameters
#' @param upper_theta Used to bound the possible values of the parameters
#' @param mcmc_iter The number of iterations of the MCMC algorithm
#' @examples

# DHM <- function(yobs, initial_theta, mcmc_iter, gibbs_cycles=10){
#
#   old_theta <- initial_theta
#
#   thetas = matrix(0,mcmc_iter,2)
#
#   for(iter in seq(1:mcmc_iter)){
#
#     # sample new parameters to propose
#     proposed_theta <- trunc_norm(old_theta, lower_theta, upper_theta)
#
#     if(any(proposed_theta) <= 0){
#
#       path[iter,]=old_theta
#
#     }else{
#
#       proposed_theta <- as.vector(proposed_theta)
#
#       # sample auxiliary variable from the likelihood at theta=proposed_theta
#       proposed_u <- simulate(~edges+triangles,coef=proposed_theta,basis=net_0,control=control.simulate(MCMC.burnin=1000,MCMC.interval=1))
#
#       # compute the approximation of the log-ratio of the partition functions
#       statistics_proposed_u <- statistics_network(proposed_u)
#       log_ratio_partitions_est <- sum(old_theta*statistics_proposed_u)-sum(proposed_theta*statistics_proposed_u)
#
#       # compute the other values appearing in the log acceptance ratio
#
#       # priors
#       log_prior_theta_old <- log_prior_theta(old_theta)
#       log_prior_theta_proposed <- log_prior_theta(proposed_theta)
#
#       # transitions  (assume symmatric for now)
#       #log_theta_old <- log_theta_kernel(proposed_theta, old_theta)
#       #log_theta_proposed <- log_theta_kernel(old_theta, proposed_theta)
#       log_theta_old <- 0
#       log_theta_proposed <- 0
#
#
#       # unnormalised likelihoods
#       yobs_stats <- statistics_network(yobs)
#       log_likelihood_old <- sum(yobs_stats*old_theta)
#       log_likelihood_proposed <- sum(yobs_stats*proposed_theta)
#
#
#       # form the acceptance rate
#       a <- log_likelihood_proposed - log_likelihood_old + log_ratio_partitions_est + log_prior_theta_proposed - log_prior_theta_old + log_theta_proposed - log_theta_old
#
#       # accept-reject
#       tmp = runif(1)
#
#       if( log(tmp) < a){
#         old_theta <- proposed_theta
#       }
#
#
#       path[iter,]=old_theta
#       if (iter%%100==0){
#         print(iter)
#         remove = c(iter:length(path[,1]))
#         plot(path[-remove,1],type="l")
#         plot(path[-remove,2],type="l")
#       }
#     }
#
#   } # end iterations
#   return(list(parameters=old_theta,estimates=path))
# }


DMH = function(Y, COV, thetas, mcmc_iter, gibbs_cycles,a){

  # Initialisation
  nCOVcols = ncol(COV)                #   number of parameters
  print(ncol(Y))
  thetaprev=numeric(nCOVcols)                  # before propose in MCMC, previous parameters
  stat = Summary(Y)
  statprop = numeric(nCOVcols)


  for(iter in seq(1:mcmc_iter)){

    if(iter%%100==0) print(iter)


    for(i in seq(1:nCOVcols)){
      thetaprev[i] = thetas[iter,i]
    }

    thetaprop = rmvnorm(n=1, mean=thetaprev, sigma=COV, method="chol")

    #proposed auxiliary variable
    U = Gibbs(Y, thetaprop, gibbs_cycles)
    statprop = Summary(U)

    #log probability ratio to determine acceptance of Outer MCMC
    #dummy = (thetaprev - thetaprop)*(statprop - stat)

    logprob = log_prior_theta(thetaprop,a)-log_prior_theta(thetaprev,a) + sum(thetaprev*statprop) - sum(thetaprev*stat) + sum(thetaprop*stat) - sum(thetaprop*statprop)

    #logprob = dummy[1]
    u = log( runif(1) )

    if( u < logprob ){
      thetas[iter+1,]=thetaprop
    }else{
      thetas[iter+1,]=thetaprev
    }

  }
  #if(iter%%100==0) print(iter)
  return(thetas)
}

Gibbs <- function(Y, theta, m){

  nrow = nrow(Y)
  star = apply(Y,1,sum)   # check what it should do...and change it
  changestat = numeric(4)

  for(l in seq(1:m)){
    for(i in seq(2:nrow)){
      for(j in seq(1:i)){ # only working on the lower triangle

        if(Y[i,j]==0){
          changestat[1] = ( choose(star[i]+1,1) + choose(star[j]+1,1) - choose(star[i],1) - choose(star[j],1) )/2 # edges   (+1 or -1)
          changestat[2] =  choose(star[i]+1,2) + choose(star[j]+1,2) - choose(star[i],2) - choose(star[j],2)      # 2-star
          changestat[3] =  choose(star[i]+1,3) + choose(star[j]+1,3) - choose(star[i],3) - choose(star[j],3)      # 3-star
          for(k in seq(1:nrow)){
            changestat[4] = changestat[4] + Y[k,i]*(Y[i,j]+1)*Y[j,k] }                                          # triangles
        }else{
          changestat[1] = ( choose(star[i],1) + choose(star[j],1) - choose(star[i]-1,1) - choose(star[j]-1,1) )/2   # edges
          changestat[2] =  choose(star[i],2) + choose(star[j],2) - choose(star[i]-1,2) - choose(star[j]-1,2)        # 2-star
          changestat[3] =  choose(star[i],3) + choose(star[j],3) - choose(star[i]-1,3) - choose(star[j]-1,3)        # 3-star
          for(k in seq(1:nrow)){
            changestat[4] = changestat[4] + Y[k,i]*Y[i,j]*Y[j,k] }                                            # triangles
        }

        r =  exp(theta*changestat)
        changestat = numeric(4)
        p = r[1]/(1+r[1])
        if( runif(1) < p  ){   # with probability  exp(theta*changestat)/(1+exp(theta*changestat))

          if(Y[i,j]==1){ star = star
          }else{
            star[i] = star[i] + 1   # reflect that we have added an edge (i->j)
            star[j] = star[j] + 1   # reflect that we have added an edge (j->i)
          }

          Y[i,j] = 1
          Y[j,i] = 1       # add the edges in the adjacency matrix
        }else{    # with probability  1   - [  exp(theta*changestat)/(1+exp(theta*changestat))   ]

          if(Y[i,j]==0){ star = star
          }else{
            star[i] = star[i] - 1   # reflect that we have removed an edge (i->j)
            star[j] = star[j] - 1    # reflect that we have removed an edge (j->i)
          }

          Y[i,j] = 0
          Y[j,i] = 0     # remove the edges in the adjacency matrix
        }

      }
    }
  }

  return(Y)
}

log_prior_theta <- function(theta,a){
  theta_1 <- dunif(theta[1],-a,a,log=TRUE)
  theta_2 <- dunif(theta[2],-a,a,log=TRUE)
  theta_3 <- dunif(theta[3],-a,a,log=TRUE)
  theta_4 <- dunif(theta[4],-a,a,log=TRUE)
  return(theta_1+theta_2+theta_3+theta_4)
}

Summary=function(A){
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
