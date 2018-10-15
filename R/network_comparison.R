#' Network comparison
#'
#' This script plots the algorithm's guess of the underlying network and compares it to the original network A
#'

#if(!require(igraph)){
 # install.packages("igraph")
#}
analyse_results = function(theta=0.5, n=50, k=10, p=0.1, alpha=0.6, beta=0.009){

  simulate <- simulate_infer(n, k, p, alpha, beta)
  g <- simulate[[1]]
  A <- simulate[[2]]
  E <- simulate[[3]]
  sim1 <- simulate[[4]][[1]]
  sim1 <- as.matrix(sim1)

  # Remove weights from simulated A
  sim1[sim1<theta]=0
  sim1[sim1>theta]=1

  # Obtain graph from the matrix provided by the algorithm
  graph <- graph_from_adjacency_matrix(sim1, mode ="undirected", weighted = NULL, diag = TRUE,
                                       add.colnames = NULL, add.rownames = NA)

  # Compare both plots
  par(mfrow=c(1,2))
  l=layout.auto(g)
  degg <- degree(g, mode="total")
  plot(g, layout=l, vertex.size=degg*2, vertex.label=NA, edge.arrow.size=0.2,vertex.color="light blue")
  deggraph <- degree(graph, mode="total")
  plot(graph, layout=l, vertex.size=deggraph*2, vertex.label=NA, edge.arrow.size=0.2,vertex.color="coral")
  par(mfrow=c(1,1))

  # Calculate accuracy metrics

  Q=sim1

  SA=A[lower.tri(A)]
  SQ=Q[lower.tri(Q)]

  FP <- sum((SQ==1) & (SA==0))
  FN <- sum((SQ==0) & (SA==1))
  TP <- sum((SQ==1) & (SA==1))
  TN <- sum((SQ==0) & (SA==0))

  FP+FN+TP+TN

  Precision <- TP / (TP+FP) #how many identified as ones are actually ones
  Recall <- TP / (TP+FN)  # proportion of ones we identify
  Accuracy <- (TP + TN) / (TP + TN + FP + FN)
  Error <- 1-Accuracy
  F_measure <- 2/ ((1 / Precision) + (1 / Recall) )

  out <- list("Precision"=Precision,"Recall"=Recall,"Accuracy"=Accuracy,"Error"=Error,"F_measure"=F_measure)
  }



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
  legend("bottomright",legend=c("F Measure", "Accuracy", "Precision", "Recall"),lty=1:4,col=c("coral", "dark green", "blue", "dark violet"), cex=1, y.intersp=2, lwd=2)
}

F_measure_plot(20, n=100)



# Plot of F_measure versus network size

F_measure_networksize <- function(k=5, n=100, alpha=0.6,beta=0.009,rho=0.1){
  size = seq(5,n,by=5)
  s <- length(size)
  Precision_vector <- 1:s
  Recall_vector <- 1:s
  F_measure_vector <- 1:s
  Accuracy_vector <- 1:s
  for (i in 1:s){
    A = sampleErdosRenyi(size[i], rho)[[1]]
    SA = A[lower.tri(A)]
    E <- interact(A,alpha,beta, size[i], k)
    simula1 <- EM(alpha0=0.4, beta0=0.02, rho0=0.15, size[i], k, E)
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
  plot(size, F_measure_vector, type="l", xlab="Network size", col="coral", ylab="",xlim=c(5,n), ylim=c(0,1),lwd=2)
  lines(size, Accuracy_vector, type="l", xlab="Network size",col="dark green", ylab="", lty=2,xlim=c(5,n) ,ylim=c(0,1),lwd=2)
  lines(size, Precision_vector, type="l", xlab="Network size",col="blue", ylab="", lty=3,xlim=c(5,n) ,ylim=c(0,1),lwd=2)
  lines(size, Recall_vector, type="l", xlab="Network size",col="dark violet", ylab="", lty=4,xlim=c(5,n) ,ylim=c(0,1),lwd=2)
  legend("bottomright",legend=c("F Measure", "Accuracy", "Precision", "Recall"),lty=1:4,col=c("coral", "dark green", "blue", "dark violet"), cex=1, y.intersp=2, lwd=2)
}


par(mfrow=c(1,1))
F_measure_plot(20, n=100)
F_measure_networksize(n=100)
