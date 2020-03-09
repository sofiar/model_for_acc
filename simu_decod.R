#### Vamos a simular un HMM autorregresivo para utilizar decoding a partir de 
#### las recurrencias de las prob forward y backward
library(ggplot2)
library(dplyr)
library(mvtnorm)
#1. Simulacion 
# 5 comportamientos. supongo todos las variables AR(1)

M=2 # Number of possible states
n=500
nrep=50

# Parameters of AR models
alphas=list()
alphas[[1]]=c(0.2,0.3,0.4) # alpha parameter for state 1 for the three acc axis
alphas[[2]]=c(0.7,0.2,-0.2) # alpha parameter for state 2 for the three acc axis
alphas[[3]]=c(-.1,-.04,0.2) # alpha parameter for state 3 for the three acc axis
alphas[[4]]=c(0.6,0.6,0.5) # alpha parameter for state 4 for the three acc axis
alphas[[5]]=c(0.4,-0.3,0.3) # alpha parameter for state 5 for the three acc axis

betas1=list() 
betas1[[1]]=c(0.4,0.5,.1) # beta1 for state 1 for the three acc axis
betas1[[2]]=c(0.2,0.35,.6) # beta1 for state 2 for the three acc axis
betas1[[3]]=c(0.5,0.6,.2) # beta1 for state 3 for the three acc axis
betas1[[4]]=c(0.1,0.2,.25) # beta1 for state 4 for the three acc axis
betas1[[5]]=c(0.3,0.4,.15) # beta1 for state 5 for the three acc axis

sigma=cbind(c(0.5,0.8,0.7),c(0.5,0.8,0.7),c(0.5,0.8,0.7),c(0.5,0.8,0.7),c(0.5,0.8,0.7))

# Initial distributions
delta=rep(1,M)/M

# para dos estados

tpm <- matrix(c(0.9, 0.1,
                0.05,0.95), byrow = T, nrow = 2)


## para 5 estados
# M=5
# tpm <- matrix(c(0.9, 0.05, 0.03, 0.01,0.01,
#                 0.1, 0.4,0.4, 0.07, 0.03,
#                 0.01, 0.15, 0.7,0.09,0.05,
#                 0.01,0.03 ,0.01, 0.61, 0.34,
#                 0.01,0.03 ,0.31, 0.11, 0.54), byrow = T, nrow = 5)

## secuencia de Estados

states <- matrix(NA, nrow = n, ncol = nrep)
for (i in 1:nrep) {
  for (j in 1:n) {
    if (j == 1) {
      states[1, i] <- sample(x = 1:M, size = 1, prob = delta)
    } else {
      states[j, i] <- sample(x = 1:M, size = 1, prob = tpm[states[j - 1, i], ])
    }
  }
}
rm(i,j)

## simu observaciones
obs <- list()
for (k in 1:nrep) obs[[k]] <- matrix(NA, nrow = n, ncol = 3)
rm(k)
# set first obs
for (k in 1:nrep){
obs[[k]][1,]=rmvnorm(1, mean = rep(0, 3), sigma = diag(sigma[,states[1,k]]))

# c(rnorm(1,0,sd=sigma[1,states[1,k]]),
#                 rnorm(1,0,sd=sigma[2,states[1,k]]),
#                 rnorm(1,0,sd=sigma[3,states[1,k]]))
}
rm(k)
for (i in 1:nrep) {
  for (j in 2:n) {
    
    mu=alphas[[states[j,i]]]+obs[[i]][(j-1),]*betas1[[states[j,i]]]
    sig=sigma[,states[j,i]]
    
    obs[[i]][j, ] = rmvnorm(n=1, mean=mu, sigma=diag(sig))
    
  }
  
}
rm(i,j)

## tidy data
tracks <- paste("T", 1:50, sep = "")
data <- data.frame(obs[[1]], state = states[, 1], timeseries = rep(tracks[1],
                                                                   500))
for (j in 2:50) data <- rbind(data, data.frame(obs[[j]], state = states[, j],
                                               timeseries = rep(tracks[j], 500)))
data$state=as.factor(data$state)
rm(j)
### plot ts 1

ggplot(data%>%filter(timeseries=='T1'),aes(x=seq(1,500),state,color=state))+geom_point()+xlab('Time')

ggplot(data%>%filter(timeseries=='T1'))+geom_line(aes(x=seq(1,500),X1),col='red')+
  geom_line(aes(x=seq(1,500),X2))+geom_line(aes(x=seq(1,500),X3),col='blue')+  
  xlab('Time')+ylab('Xs')

######################################################################
####### Generamos las entradas para el decoding ######################
######################################################################

# M = number of states
# n = lenght of the time series
# log_gamma = transition probability matrix for a homogeneous state process
#             log transformed 
#             :: MxM matrix
# log_allprobs = logarithm of the state-dependent densities evaluated 
#                at each observation 
#                :: n x M matrix
# log_delta  = log initial distribution 
#              :: vector of length m



log_allprobs=list()

for (i in 1:nrep)
  {log_allprobs[[i]]=matrix(NA,ncol=M,nrow=n)}

rm(i)
for (i in 1:nrep) {
  for (m in 1:M) {
    # prob at time 1
    log_allprobs[[i]][1,m]=dmvnorm(obs[[i]][1,], mean=c(0,0,0), sigma=diag(sigma[,m]), log = TRUE)
        for (t in 2:n)  
        {
        mu=alphas[[m]]+obs[[i]][t-1,]*betas1[[m]]
        sig=sigma[,m]
        log_allprobs[[i]][t,m]= dmvnorm(obs[[i]][t,], mean=mu, sigma=diag(sig), log = TRUE)
        }
}
}
rm(i,m,t)



log_delta=log(delta)
log_gamma=log(tpm)



