#' Network comparison
#'
#' This script plots the algorithm's guess of the underlying network and compares it to the original network A
#'

if(!require(igraph)){
  install.packages("somepackage")
}

g <- sampleErdosRenyi(n=50,p=0.1)[[2]]
sim1 <- as.matrix(simulation1)[[1]]

# Remove weights from simulated A

sim1[sim1<0.3]=0
sim1[sim1>0.8]=1

# Check that the entire matrix has been converted to 1s and 0s

sum((sim1==0)|(sim1==1))

# Obtain graph from the matrix provided by the algorithm

graph <- graph_from_adjacency_matrix(sim1, mode ="undirected", weighted = NULL, diag = TRUE,
                                     add.colnames = NULL, add.rownames = NA)

# Compare both plots

par(mfrow=c(1,2))
l=layout=layout.auto(g)
degg <- degree(g, mode="total")
plot(g, layout=l, vertex.size=degg*2, vertex.label=NA, edge.arrow.size=0.2,vertex.color="light blue")
deggraph <- degree(graph, mode="total")
plot(graph, layout=l, vertex.size=deggraph*2, vertex.label=NA, edge.arrow.size=0.2,vertex.color="coral")
par(mfrow=c(1,1))

# Calculate proportion of edges correctly identified by the algorithm

proportion_guess <- (sum(A==sim1))/prod(dim(A)) #not very useful cause most entries are zero so easy to guess
proportion_ofedges_guess <- sum(A[A==1]==sim1[A==1])/sum(A)

proportion_guess
proportion_ofedges_guess


