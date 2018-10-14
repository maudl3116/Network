#' Network comparison
#'
#' This script plots the algorithm's guess of the underlying network and compares it to the original network A
#'

#if(!require(igraph)){
 # install.packages("igraph")
#}


# Plot of F_measure versus observed days k

F_measure_plot <- function(k, n, alpha=0.6,beta=0.009,rho=0.1){
  A = sampleErdosRenyi(n, rho)[[1]]
  SA = A[lower.tri(A)]
  Precision_vector <- 1:k
  Recall_vector <- 1:k
  F_measure_vector <- 1:k
  Accuracy_vector <- 1:k
  for (i in 1:k){
    E <- interact(A,alpha,beta, n,i)
    simula1 <- EM(alpha0=0.4, beta0=0.02, rho0=0.15, n, i, E)
    Q <- as.matrix(simula1)[[1]]
    SQ=Q[lower.tri(Q)]

    theta = 0.5
    SQ[SQ<theta]=0
    SQ[SQ>theta]=1

    FP <- sum((SQ==1) & (SA==0))
    FN <- sum((SQ==0) & (SA==1))
    TP <- sum((SQ==1) & (SA==1))
    TN <- sum((SQ==0) & (SA==0))
    Accuracy_vector[i] <- (TP + TN) / (TP + TN + FP + FN)
    Precision_vector[i] <- TP / (TP+FP)
    Recall_vector[i] <- TP / (TP+FN)
    F_measure_vector[i] <- 2/ ((1 / Precision_vector[i]) + (1 / Recall_vector[i]))
  }
  plot(F_measure_vector, type="l", xlab="Number of days", col="coral", ylab="",xlim=c(1,k), ylim=c(0,1),lwd=2)
  lines(Accuracy_vector, type="l", xlab="Number of days",col="dark green", ylab="", lty=2,xlim=c(1,k) ,ylim=c(0,1),lwd=2)
  lines(Precision_vector, type="l", xlab="Number of days",col="blue", ylab="", lty=3,xlim=c(1,k) ,ylim=c(0,1),lwd=2)
  lines(Recall_vector, type="l", xlab="Number of days",col="dark violet", ylab="", lty=4,xlim=c(1,k) ,ylim=c(0,1),lwd=2)
  legend("bottomright",legend=c("F Measure", "Accuracy", "Precision", "Recall"),lty=1:4,col=c("coral","dark green", "blue", "dark violet"), cex=1, y.intersp=0.6, lwd=2)
}

F_measure_plot(20, n=100)

