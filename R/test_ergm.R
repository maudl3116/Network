data(florentine)
plot(flobusiness)

# Parameter estimation via MLE or MPLE #
formula <- flobusiness ~ edges + kstar(2) + kstar(3) + triangle
m <-ergm(formula,estimate="MPLE")

# summary statistics for data X
X = as.matrix(flobusiness)

# initialisation
initial= m$coef
sigma = 0.15
mcmc_iter = 10000
thetas = matrix(0,mcmc_iter,4)
thetas[1,]=runif(4,-5,5)

ptm <- proc.time()
Result = exchange(X, diag(sigma^2,4,4), mcmc_iter, thetas, 5)

ptme <- proc.time()
ptme - ptm

# analysis results

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


#Resultsummary = rbind( t(bmmat(Result)), HPDinterval(as.mcmc(Result), prob = 0.95)[,1],
#                     HPDinterval(as.mcmc(Result), prob = 0.95)[,2], c(ess(Result[,1]),ess(Result[,2]),ess(Result[,3]),ess(Result[,4])),
#                    rep(length(unique(Result[,1]))/outer,4), rep(Ltime,4))


Resultsummary = rbind( t(bmmat(Result)), c(ess(Result[,1]),ess(Result[,2]),ess(Result[,3]),ess(Result[,4])),
                      rep(length(unique(Result[,1]))/outer,4), rep(Ltime,4))

save(Resultsummary ,Result, file = "ExchangeErgm_10000_2.RData")


