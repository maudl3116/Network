
install.packages("igraph")
library(igraph)

par(mfrow=c(1,4))
g <- erdos.renyi.game(2, 1, type = c("gnp"), directed=FALSE, loops=FALSE)
plot(g, vertex.size=15, vertex.label=NA, edge.arrow.size=20,vertex.color="light blue")
g <- erdos.renyi.game(3, 0.6, type = c("gnp"), directed=FALSE, loops=FALSE)
plot(g, vertex.size=15, vertex.label=NA, edge.arrow.size=20,vertex.color="light blue")
g <- erdos.renyi.game(4, 0.4, type = c("gnp"), directed=FALSE, loops=FALSE)
plot(g, vertex.size=15, vertex.label=NA, edge.arrow.size=0.2,vertex.color="light blue")
g <- erdos.renyi.game(3, 1, type = c("gnp"), directed=FALSE, loops=FALSE)
plot(g, vertex.size=15, vertex.label=NA, edge.arrow.size=0.2,vertex.color="light blue")
par(mfrow=c(1,1))
