# generating an observation yobs from the ground truth network
n=10
m <- matrix(rbinom(n*n,1,.5),n,n)
diag(m) <- 0
net_0=network(m,directed=FALSE)
yobs = simulate(~edges+triangles,coef=c(3,5),basis=net_0,control=control.simulate(MCMC.burnin=1000,MCMC.interval=1))


out = ergm_no_noise_fit(n, yobs, c(1,7), c(-Inf, -Inf), c(Inf, Inf), 10000)
burn = c(seq(1:900))
plot(out$estimates[-burn,1],type="l")
plot(out$estimates[-burn,2],type="l")
hist(out$estimates[-burn,1],breaks=20)
hist(out$estimates[-burn,2],breaks=20)
