#' Network comparison
#'
#' This script plots the algorithm's guess of the underlying network and compares it to the original network A
#'

if(!require(igraph)){
  install.packages("igraph")
}

#g <- sampleErdosRenyi(n=50,p=0.1)[[2]]
sim1 <- as.matrix(simulation1)[[1]]

# Remove weights from simulated A
theta = 0.5
sim1[sim1<theta]=0
sim1[sim1>theta]=1

# Check that the entire matrix has been converted to 1s and 0s

sum((sim1==0)|(sim1==1))

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

# Calculate proportion of edges correctly identified by the algorithm

Q=sim1

SA=A[lower.tri(A)]
SQ=Q[lower.tri(Q)]

FP <- sum((SQ==1) & (SA==0))
FN <- sum((SQ==0) & (SA==1))
TP <- sum((SQ==1) & (SA==1))
TN <- sum((SQ==0) & (SA==0))

FP+FN+TP+TN



