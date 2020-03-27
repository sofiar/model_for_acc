
x = sim$obs
od = 'norm'
rd = 'pois'



## Aca arranca la funcion hsmm.smoth
  
inputData <- as.vector(x)
tau       <- length(x) 
error     <- as.integer(0)

# determine number of states
J <- length(delta)

# selection of the maximum runlength
M = as.integer(min(tau, 1000)) # este 1000 es arbitrario


# variables calculated by Forward-Backward alg.
F    <- as.double(rep(0, times = J * tau))
L    <- as.double(rep(0, times = J * tau))
G    <- as.double(rep(0, times = J * tau))
L1   <- as.double(rep(0, times = J * tau))
N    <- as.double(rep(0, times = tau))
Norm <- as.double(rep(0, times = J * tau))
eta  <- as.double(rep(0, times = J * M))
xi   <- as.double(rep(0, times = J * M))

# Store variables for calling FB
FB.p.tpm      <- t(tpm)
dim(FB.p.tpm) <- c(J * J)
FB.pi.ini     <- delta

FB.d          <- get.d(rd, J, M, param = lambdas)
FB.pdf        <- get.pdf(inputData, od, J, M, param = odpar)   

# Call Forward-Backward Alg.
FB.result <- FB(0, tau, J, M, FB.d, 
                FB.p.tpm, FB.pi.ini, FB.pdf, F, L, G, 
                L1, N, Norm, eta, xi, error)


L=FB.result[[10]]
F= FB.result[[9]]
G=FB.result[[11]]
N=FB.result[[13]]
L1= FB.result[[12]]
# calculate llh 
llh = sum(log(N[1:tau]))


# pasa a matriz
dim(L)= c(tau, J)
dim(F)= c(tau, J)
dim(G)= c(tau, J)
dim(L1)= c(tau, J)

## 0-1 LF
pred2=apply(L, 1, which.max)
Error.pred=sum(pred2==sim$path)/length(x)
Error.pred


## No entiendo como hace para calcular el valor de L
## pruebo a ver si es el producto de las probas


L.sofi=matrix(NA,ncol=J,nrow=tau)
L.sofi[tau,]=F[tau,]



for (t in (tau-1):1)
{
  for (q in 1:J)
  {
    L.sofi[t,q]=L.sofi[t+1,q]+
      F[t,q]*sum(G[t+1,]*tpm[q,])- # porque la t(tpm)?
      G[t+1,q]*sum(F[t,]*tpm[,q])
  }
}

sum(L.sofi!=L)
# SOn iguales 
pred.sofi=apply(L.sofi, 1, which.max)

sum(pred2==pred.sofi)/length(x)

sum(pred.sofi==sim$path)/length(x)


######
library(Rcpp)
Rcpp::sourceCpp('prueba.cpp')
prueba(d[1,],StateIn[j,],pdf[j,],pi[j],J,N)




F[j,t]=0
estr=sum(F[,t-1]*p[,j])
pobs=pi[j]*pdf[j,1]*pdf[j,2]/(N[2]*N[1])
F[j,t]=F[j,t]+estr*d[j,i]*pobs  
F[j,t]
 
