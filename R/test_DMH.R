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

#### Conduct Liang's DMH ###
outer = 10000  # was 30,000
cycle = 20                   # inner sampler is cycle*N(N-1)/2

COV = diag(0.01,4)
a=5
initial = matrix(runif(4,-a,a),outer+1,4)


ptm <- proc.time()
Liang = DMH(X, COV, initial, outer, cycle,a)
ptme <- proc.time()
ptme - ptm


# analysis results

# qualitative analysis
par(mfrow=c(2,2))
ts.plot(Liang[,1])
ts.plot(Liang[,2])
ts.plot(Liang[,3])
ts.plot(Liang[,4])

par(mfrow=c(2,2))
hist(Liang[,1])
hist(Liang[,2])
hist(Liang[,3])
hist(Liang[,4])

# quantitative analaysis
bmmat(Liang)  # applies batch means for each parameter trace
#HPDinterval(as.mcmc(Liang), prob = 0.95)   # ? problem highest probability density intervals
length(unique(Liang[,1]))/outer
ess(Liang[,1])  # effective sample size
Ltime = (ptme - ptm)[[1]]


# Liangsummary = rbind( t(bmmat(Liang)), HPDinterval(as.mcmc(Liang), prob = 0.95)[,1],
#                       HPDinterval(as.mcmc(Liang), prob = 0.95)[,2], c(ess(Liang[,1]),ess(Liang[,2]),ess(Liang[,3]),ess(Liang[,4])),
#                       rep(length(unique(Liang[,1]))/outer,4), rep(Ltime,4))
Liangsummary = rbind( t(bmmat(Liang)), c(ess(Liang[,1]),ess(Liang[,2]),ess(Liang[,3]),ess(Liang[,4])),
                      rep(length(unique(Liang[,1]))/outer,4), rep(Ltime,4))
save(Liangsummary ,Liang, file = "DMH_10000.RData")


bmmat(Liang)


