#practice SEM modeling

library("OpenMx", lib.loc="~/R/win-library/3.1")

#set.seed(555)

load  <- matrix( c(1,0, 1,0, 1,0, 0,1, 0,1, 0,1), 
                 byrow=TRUE, nrow=6, ncol=2)
psi   <- matrix(c(4, 2, 2, 3), 2, 2)
theta <- diag(c(.2, .2, .2, .2, .2, .2))
cov.hat <- load %*% psi %*% t(load) + theta

L.hat <- chol(cov.hat)

Y <- matrix(rnorm(6000),1000,6) %*% t(L.hat)

colnames(Y) <- c("y1","y2","y3","y4","y5","y6")
Y.cov <- cov(Y)

manifestVars <- c("y1","y2","y3","y4","y5","y6")
latentVars <- c("F1", "F2")

covdata <- mxData(observed=Y.cov, type="cov", numObs=1000)

a.mat <- mxMatrix(type="Full", nrow=8, ncol=8, byrow=TRUE, name="A", 
                  free = c(F,F,F,F,F,F,T,F,
                           F,F,F,F,F,F,T,F,
                           F,F,F,F,F,F,T,F,
                           F,F,F,F,F,F,F,T,
                           F,F,F,F,F,F,F,T,
                           F,F,F,F,F,F,F,T,
                           F,F,F,F,F,F,F,F,
                           F,F,F,F,F,F,F,F), 
                  values = c(0,0,0,0,0,0,1,0,
                             0,0,0,0,0,0,1,0,
                             0,0,0,0,0,0,1,0,
                             0,0,0,0,0,0,0,1,
                             0,0,0,0,0,0,0,1,
                             0,0,0,0,0,0,0,1,
                             0,0,0,0,0,0,0,0,
                             0,0,0,0,0,0,0,0),
                  labels = c(NA,NA,NA,NA,NA,NA,"L11", NA,
                             NA,NA,NA,NA,NA,NA,"L21", NA,
                             NA,NA,NA,NA,NA,NA,"L31", NA,
                             NA,NA,NA,NA,NA,NA, NA,  "L42",
                             NA,NA,NA,NA,NA,NA, NA,  "L52",
                             NA,NA,NA,NA,NA,NA, NA,  "L62",
                             NA,NA,NA,NA,NA,NA, NA,   NA,
                             NA,NA,NA,NA,NA,NA, NA,   NA))

s.mat <- mxMatrix(type="Symm", nrow=8, ncol=8, byrow=TRUE, name="S",
                  free = c(T,F,F,F,F,F,F,F, 
                           F,T,F,F,F,F,F,F, 
                           F,F,T,F,F,F,F,F, 
                           F,F,F,T,F,F,F,F, 
                           F,F,F,F,T,F,F,F, 
                           F,F,F,F,F,T,F,F, 
                           F,F,F,F,F,F,F,T, 
                           F,F,F,F,F,F,T,F),     
                  values = c(.2, 0, 0, 0, 0, 0, 0, 0,
                             0,.2, 0, 0, 0, 0, 0, 0,
                             0, 0,.2, 0, 0, 0, 0, 0,
                             0, 0, 0,.2, 0, 0, 0, 0,
                             0, 0, 0, 0,.2, 0, 0, 0,
                             0, 0, 0, 0, 0,.2, 0, 0,
                             0, 0, 0, 0, 0, 0, 4, 2,
                             0, 0, 0, 0, 0, 0, 2, 3), 
                  labels = c("e1",NA,NA,NA,NA,NA,NA,NA, 
                             NA,"e2",NA,NA,NA,NA,NA,NA,
                             NA,NA,"e3",NA,NA,NA,NA,NA,
                             NA,NA,NA,"e4",NA,NA,NA,NA,
                             NA,NA,NA,NA,"e5",NA,NA,NA,
                             NA,NA,NA,NA,NA,"e6",NA,NA,
                             NA,NA,NA,NA,NA,NA,"P11","P21",
                             NA,NA,NA,NA,NA,NA,"P21","P22"))

f.mat <- mxMatrix(type="Full", nrow=6, ncol=8, free=FALSE, byrow=TRUE, name="F",
                  values = c(1,0,0,0,0,0,0,0, 
                             0,1,0,0,0,0,0,0,
                             0,0,1,0,0,0,0,0,
                             0,0,0,1,0,0,0,0,
                             0,0,0,0,1,0,0,0,
                             0,0,0,0,0,1,0,0))

exp <- mxExpectationRAM("A", "S", "F", 
                        dimnames=c(manifestVars, latentVars)) 

funML <- mxFitFunctionML()

diverto <- mxModel("Let us have some fun with R", 
                   covdata, a.mat, s.mat, f.mat, exp, funML)

diverto.run <- mxRun(diverto)     