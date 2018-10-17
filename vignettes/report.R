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

