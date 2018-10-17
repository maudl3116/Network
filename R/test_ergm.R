# # generating an observation yobs from the ground truth network
# n=10
# m <- matrix(rbinom(n*n,1,.5),n,n)
# diag(m) <- 0
# net_0=network(m,directed=FALSE)
# yobs = simulate(~edges+triangles,coef=c(3,5),basis=net_0,control=control.simulate(MCMC.burnin=1000,MCMC.interval=1))
#
#
# out = ergm_no_noise_fit(n, yobs, c(1,7), c(-Inf, -Inf), c(Inf, Inf), 10000)
# burn = c(seq(1:900))
# plot(out$estimates[-burn,1],type="l")
# plot(out$estimates[-burn,2],type="l")
# hist(out$estimates[-burn,1],breaks=20)
# hist(out$estimates[-burn,2],breaks=20)



data(florentine)
plot(flobusiness)

# Parameter estimation via MLE or MPLE #
formula <- flobusiness ~ edges + kstar(2) + kstar(3) + triangle
m <-ergm(formula,estimate="MPLE")



# summary statistics for data X
X = as.matrix(flobusiness)


# initialisation
initial= m$coef
sigma = 0.5
mcmc_iter = 100
thetas = matrix(0,mcmc_iter,4)
thetas[1,]=initial

ptm <- proc.time()
Liang = exchange(X, diag(sigma^2,4,4), mcmc_iter, thetas, 5)

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
HPDinterval(as.mcmc(Liang), prob = 0.95)   # ? problem highest probability density intervals
length(unique(Liang[,1]))/outer
ess(Liang[,1])  # effective sample size
Ltime = (ptme - ptm)[[1]]


Liangsummary = rbind( t(bmmat(Liang)), HPDinterval(as.mcmc(Liang), prob = 0.95)[,1],
                      HPDinterval(as.mcmc(Liang), prob = 0.95)[,2], c(ess(Liang[,1]),ess(Liang[,2]),ess(Liang[,3]),ess(Liang[,4])),
                      rep(length(unique(Liang[,1]))/outer,4), rep(Ltime,4))

save(Liangsummary ,Liang, file = "ExchangeErgm.RData")



