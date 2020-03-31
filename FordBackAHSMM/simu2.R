# Simulacion a mano sin usar el paquete 

library(ggplot2)
library(dplyr)
library(mvtnorm)
library(glmnet)
library(hsmm)
source('simFun.R')
source('FBA.R')
source('FBcpp.R')

### 1.Parametros de las simulaciones  
delta  <- rep(1/3, 3)
tpm <- t(matrix(c(0, 0.5, 0.5,
                0.7, 0, 0.3,
                0.8, 0.2, 0), 3, byrow = TRUE))
lambdas  <- c(3,6,8)
odpar  <- list(mean = c(-1.5, 0, 1.5), var = c(0.5, 0.6, 0.8))
n = 200

# simulation and calculation of the prediction error
Nrep=1000
E.pred=numeric(Nrep)
seeds=sample(seq(1,1e7),size=Nrep,replace=F)
for (l in 1:Nrep)
{
sim=sim.ts.NP(delta, lambdas, tpm, params=list(mu= c(-1.5, 0, 1.5),
                                           sd=sqrt(c(0.5, 0.6, 0.8))),
                                           n, seed=seeds[l],M=1000)
obs=sim$Obs
States=sim$States

#################################################################
#################         4. Prediccion       ###################
#################################################################
### Parametros para la prediccion
###
### inputData : vector de observaciones
### tau: largo de la serie
### J: numero de estados posibles
### M: selection of the maximum runlength
### FB.d: probs of sojuorn times
### FB.pdf: vector of allprobs (pdfs)
### tpm.predict: tpm as vector

J = length(delta)
tau = length(obs) 
M = as.integer(min(tau, 1000))

## Sojourn time 
FB.d=matrix(NA,M,J)
for (j in 1:J)
{
  FB.d[,j]= c(dpois(c(0:(M - 1)), lambdas[j]))
}
FB.d=t(FB.d)
## PDF
FB.pdf=matrix(NA,J,tau)
for (j in 1:J)
{
FB.pdf[j,] <- dnorm(obs[1:tau], mean = odpar$mean[j], sd = sqrt(odpar$var[j]))
}
lower_bound <- 1e-300
FB.pdf[FB.pdf < lower_bound] <- lower_bound

FB.Result2=FBAcpp(pi=delta,pdf=FB.pdf,d=FB.d,p=tpm)
if (abs(sum(FB.Result2$Forwrd[,tau])-1)>0.00001)
{
  break
  print(sum(FB.Result2$Forwrd[,tau]))
}

### Prediccion
pred2=apply(FB.Result2$Gamma,2,which.max)
E.pred[l]=mean(pred2==States)
print(l)
}

hist(E.pred)
