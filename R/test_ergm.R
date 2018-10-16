n=10
m <- matrix(rbinom(n*n,1,.5),n,n)
diag(m) <- 0
net_0=network(m,directed=FALSE)
yobs = simulate(~edges+triangles,coef=c(-3,3),basis=net_0,control=control.simulate(MCMC.burnin=100,MCMC.interval=1000))


out = ergm_no_noise_fit(n, yobs, c(-1,10), c(-Inf, -Inf), c(Inf, Inf), 10000)
plot(out$estimates[,1],type="l")
plot(out$estimates[,2],type="l")
out
