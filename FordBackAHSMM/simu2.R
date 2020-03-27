# Simulacion a mano sin usar el paquete 

library(ggplot2)
library(dplyr)
library(mvtnorm)
library(glmnet)
library(hsmm)
source('funAuxhsmm.R')
### 1.Parametros de las simulaciones  
delta  <- rep(1/3, 3)
tpm <- t(matrix(c(0, 0.5, 0.5,
                0.7, 0, 0.3,
                0.8, 0.2, 0), 3, byrow = TRUE))
lambdas  <- c(3,6,8)
odpar  <- list(mean = c(-1.5, 0, 1.5), var = c(0.5, 0.6, 0.8))

n = 300
J = 3 # number of possible states

### 2. Simulation of the NAS sequence 
set.seed(999)
states=numeric(n)
for (j in 1:n) {
  if (j == 1) {
    states[1] <- sample(x = 1:J, size = 1, prob = delta)
  } else {
    states[j] <- sample(x = 1:J, size = 1, prob = tpm[states[j -1], ])
  }
}

### 3. Simultation of the sojourn times 
Stimes=numeric(n)

RL=matrix(get.d(rd='pois',J,M=1000,param=list(lambda=lambdas)),nrow=J,byrow=TRUE)

i=states[1]
for (j in 1:n) {
  #Stimes[j] <- rpois(1,lambdas [states[j]])
  Stimes[j] =sample(c(1:length(RL[i, RL[i,]!=0])), prob=RL[i, RL[i,]!=0], 1)  
}
totL=sum(Stimes)

States=rep(states,times=Stimes)
## 3. Observations (Nomales independientes)

obs <- rep(NA,totL)

# set first obs
for (j in 1:totL) {
    mu=odpar$mean[States[j]]
    sig=sqrt(odpar$var[States[j]])
    obs[j] = rnorm(n=1, mean=mu, sd=sig)
  }


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

source('FBA.R')
source('FBcpp.R')

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

FB.Result=FBA(pi=delta,pdf=FB.pdf,d=FB.d,p=tpm)
FB.Result2=FBAcpp(pi=delta,pdf=FB.pdf,d=FB.d,p=tpm)

### Prediccion
pred=apply(FB.Result$Gamma,2,which.max)
mean(pred==States)


pred2=apply(FB.Result2$Gamma,2,which.max)
mean(pred2==States)


