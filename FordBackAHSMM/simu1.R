### Voy a debuggear el codigo 

library(ggplot2)
library(dplyr)
library(mvtnorm)
library(glmnet)
library(hsmm)
source('funAuxhsmm.R')

# Simulating observations: 
delta  <- rep(1/3, 3)
tpm <- matrix(c(0, 0.5, 0.5,
                   0.7, 0, 0.3,
                   0.8, 0.2, 0), 3, byrow = TRUE)
lambdas  <- list(lambda = c(3,6,8))
odpar  <- list(mean = c(-1.5, 0, 1.5), var = c(0.5, 0.6, 0.8))

# en realidad es la t(tpm) nose porque la mete asi. comprobar con simuaciones
sim    <- hsmm.sim(n = 300, od = "norm", rd = "pois", 
                   pi.par = delta, tpm.par = t(tpm), 
                   rd.par = lambdas, od.par = odpar, seed = 1606) # En realidad en es la t(tpm)


#################################################################
#################         4. Prediccion       ###################
#################################################################


inputData = sim$obs
tau= length(inputData)
error= as.integer(0)
od = 'norm'
rd = 'pois'
J = length(delta)
M = as.integer(min(tau, 1000)) # este 1000 es arbitrario

# Store variables for calling FB
FB.p.tpm      <- t(tpm)
dim(FB.p.tpm) <- c(J * J)
FB.d          <- get.d(rd, J, M, param = lambdas)
FB.pdf        <- get.pdf(inputData, od, J, M, param = odpar)   




# variables calculated by Forward-Backward alg.
F    <- as.double(rep(0, times = J * tau))
L    <- as.double(rep(0, times = J * tau))
G    <- as.double(rep(0, times = J * tau))
L1   <- as.double(rep(0, times = J * tau))
N    <- as.double(rep(0, times = tau))
Norm <- as.double(rep(0, times = J * tau))
eta  <- as.double(rep(0, times = J * M))
xi   <- as.double(rep(0, times = J * M))


# Call Forward-Backward Alg.
FB.result <- FB(0, tau, J, M, FB.d, 
                FB.p.tpm, delta, FB.pdf, F, L, G, 
                L1, N, Norm, eta, xi, error)


Lsim1=FB.result[[10]]
Fsim1= FB.result[[9]]
Gsim1=FB.result[[11]]
Nsim1=FB.result[[13]]
L1sim1= FB.result[[12]]
Nsim1 = FB.result[[13]]
Normsim1 <- FB.result[[14]]

# pasa a matriz
dim(Lsim1)= c(tau, J)
dim(Fsim1)= c(tau, J)
dim(Gsim1)= c(tau, J)
dim(L1sim1)= c(tau, J)
dim(Normsim1)= c(tau, J)

Fsim1=t(Fsim1)
Gsim1=t(Gsim1)
Normsim1=t(Normsim1)


dim(Fsim1)= c(J,tau)

## 0-1 LF
pred=apply(Lsim1, 1, which.max)
Error.pred=sum(pred==sim$path)/length(sim$obs)
Error.pred




