data(florentine)
plot(flobusiness)

# Parameter estimation via MLE or MPLE #
formula <- flobusiness ~ edges + kstar(2) + kstar(3) + triangle
m <-ergm(formula,estimate="MPLE")
m$coef

# summary statistics for data X
X = as.matrix(flobusiness)
g=graph_from_adjacency_matrix(X)
plot(g, layout=layout.auto, vertex.size=6, vertex.label=NA, edge.arrow.size=0.2,vertex.color="light blue")
Summary(X)
summary(formula)

#### Conduct Result's DMH ###
outer = 10000  # was 30,000
cycle = 20                   # inner sampler is cycle*N(N-1)/2

COV = diag(0.01,4)
a=5
initial = matrix(runif(4,-a,a),outer+1,4)


ptm <- proc.time()
Result = DMH(X, COV, initial, outer, cycle,a)
ptme <- proc.time()
ptme - ptm


# analysis Result

# qualitative analysis
par(mfrow=c(2,2))
ts.plot(Result[,1])
ts.plot(Result[,2])
ts.plot(Result[,3])
ts.plot(Result[,4])

par(mfrow=c(2,2))
hist(Result[,1])
hist(Result[,2])
hist(Result[,3])
hist(Result[,4])

# quantitative analaysis
bmmat(Result)  # applies batch means for each parameter trace
#HPDinterval(as.mcmc(Result), prob = 0.95)   # ? problem highest probability density intervals
length(unique(Result[,1]))/outer
ess(Result[,1])  # effective sample size
Ltime = (ptme - ptm)[[1]]


# Resultsummary = rbind( t(bmmat(Result)), HPDinterval(as.mcmc(Result), prob = 0.95)[,1],
#                       HPDinterval(as.mcmc(Result), prob = 0.95)[,2], c(ess(Result[,1]),ess(Result[,2]),ess(Result[,3]),ess(Result[,4])),
#                       rep(length(unique(Result[,1]))/outer,4), rep(Ltime,4))
Resultsummary = rbind( t(bmmat(Result)), c(ess(Result[,1]),ess(Result[,2]),ess(Result[,3]),ess(Result[,4])),
                      rep(length(unique(Result[,1]))/outer,4), rep(Ltime,4))
save(Resultsummary ,Result, file = "DMH_10000_2.RData")



