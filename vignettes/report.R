## ----pressure, echo=FALSE, fig.pos="H", out.width = '65%',fig.cap="\\label{fig:graph_model}Graphical Model"----
knitr::include_graphics("GM2.png")


## ----echo=FALSE, fig.align='center', fig.pos="H", out.width = '60%',out.height='60%',fig.cap="\\label{fig:triangle}Example of edge, 2-star, 3-star and triangle"----
knitr::include_graphics("triangle.png")


## ----echo=FALSE, fig.align='center', fig.pos="H", out.width = '20%',out.height='20%',fig.cap="\\label{fig:triangle}The augmented model"----
knitr::include_graphics("exchange.png")


## ----echo=FALSE, results='hide',fig.keep='none',fig.pos="H"--------------
   # generate ground truth network
  output <- sampleErdosRenyi(n,rho)
  g <- output[[2]]
  A <- output[[1]]

  # generate noisy observations of the ground truth network
  E <- interact(A,alpha,beta, n,k)
  simulation <- EM(alpha0=0.4, beta0=0.02, rho0=0.15, n, k, E)

## ----echo=FALSE, fig.pos="H",fig.cap="\\label{fig:comparison}(left) Ground truth underlying network (right) Inferred underlying network"----
  out=analyse_results(t,n,k,rho,alpha,beta)

## ----echo=FALSE,fig.pos="H",fig.height = 2,fig.width=8,fig.align = "center",fig.cap="\\label{fig:figs} Convergence of the parameter estimates"----

m1<-simulation[[6]][c(1:simulation[[5]]),]
par(mfrow = c(1, 3))
plot(m1[,1], type="l",xlab="iteration number",ylab=expression(alpha))
plot(m1[,2],type="l",xlab="iteration number",ylab=expression(beta))
plot(m1[,3],type="l",xlab="iteration number",ylab=expression(rho))

## ----echo=FALSE,fig.pos="H",fig.height = 4,fig.width=10,fig.align = "center",fig.cap="\\label{fig:figs}Performance metrics versus the number of repeated observations k"----
par(mfrow=c(1,2))
F_measure_plot(20, n=100)
F_measure_networksize(n=100)

