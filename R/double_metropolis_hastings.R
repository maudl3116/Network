#' This function makes it possible to fit an ERGM model (i.e. sample from the posterior distribution on the parameters) using an enhanced MCMC technique: the Double Metropolis Hastings Algorithm.
#' @param Y one observation arising from the unknown network structure (must be an adjacency matrix)
#' @param COV scaling of the random walk.
#' @param thetas a matrix storing the updates of the parameters for each DMH iteration. The first line of the matrix contains the initial values of the parameters.
#' @param mcmc_iter The number of iterations of the outer Metropolis Hastings sampler
#' @param gibbs_cycles The number of iterations of the inner sampler (the number of Gibbs updates before considering the last state as a draw from the likelihood P(x|theta) for the auxiliary variable sample)
#' @param a the scale of the uniform prior for each parameter

DMH = function(Y, COV, thetas, mcmc_iter, gibbs_cycles,a){

  # Initialisation
  nCOVcols = ncol(COV)                         # gives the number of parameters in the model
  thetaprev=numeric(nCOVcols)                  # before propose in MCMC, previous parameters
  stat = Summary(Y)                            # compute the sufficient statistics for the observation
  statprop = numeric(nCOVcols)                 # initialisation of the sufficient statistics of the proposed theta


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

    logprob = log_prior_theta(thetaprop,a)-log_prior_theta(thetaprev,a) + sum(thetaprev*statprop) - sum(thetaprev*stat) + sum(thetaprop*stat) - sum(thetaprop*statprop)

    u = log( runif(1) )

    if( u < logprob ){
      thetas[iter+1,]=thetaprop
    }else{
      thetas[iter+1,]=thetaprev
    }

  }
  return(thetas)
}

#' This function performs m full Gibbs update of an initial graph Y
#' @param m The number of full pass over Y's components, i.e. the number of Gibbs updates
#' @param Y The initial graph which is updated m times with Gibbs, to get a draw from the likelihood
#' @param theta The parameters of the likelihood.
Gibbs <- function(Y, theta, m){

  nrow = nrow(Y)
  star = apply(Y,1,sum)    # sum of all the lines in Y, used to compute the sufficient statistics of Y
  changestat = numeric(4)

  for(l in seq(1:m)){
    for(i in seq(2:nrow)){
      for(j in seq(1:i)){ # only working on the lower triangle

        if(Y[i,j]==0){  # if we change Y[i,j] from 0 to 1, how will it impact te sufficient statistics ?
          changestat[1] = ( choose(star[i]+1,1) + choose(star[j]+1,1) - choose(star[i],1) - choose(star[j],1) )/2 # edges   (+1 or -1)
          changestat[2] =  choose(star[i]+1,2) + choose(star[j]+1,2) - choose(star[i],2) - choose(star[j],2)      # 2-star
          changestat[3] =  choose(star[i]+1,3) + choose(star[j]+1,3) - choose(star[i],3) - choose(star[j],3)      # 3-star
          for(k in seq(1:nrow)){
            changestat[4] = changestat[4] + Y[k,i]*(Y[i,j]+1)*Y[j,k] }                                          # triangles
        }else{   # if we change Y[i,j] from 1 to 0, how will it impact te sufficient statistics ?
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

#' This function computes the statistics of the graph A
#' @param A The adjacency matrix of the graph of which we want to compute the statistics
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

