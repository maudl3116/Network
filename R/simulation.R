#' Simulating
#'
#' This function allows you to simulate network data and perform inference via the Expectation Maximization algorithm
simulate_infer=function(n=50,k=10,p=0.1,alpha=0.6, beta=0.009){

  # generate ground truth network
  output <- sampleErdosRenyi(n,p)
  g <- output[[2]]
  A <- output[[1]]

  # generate noisy observations of the ground truth network
  E <- interact(A,alpha,beta, n,k)

  # perform inference via the EM algorithm
  simulation <- EM(alpha0=0.4, beta0=0.02, rho0=0.15, n, k, E)
}



analyse_results = function(theta=0.5){

  sim1 <- as.matrix(simulation)[[1]]

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
