########################################################################################
################################ Inference with STAN ###################################
########################################################################################

library(rstan)
library(tidyverse)

#############################
## Toy example HMM 2States ##
#############################

# Simulation

# Number of states
N <- 2
# transition probabilities
Gamma <- matrix(c(0.9,0.1,0.1,0.9),2,2)
# initial distribution set to the stationary distribution
delta <- solve(t(diag(N)-Gamma +1), rep(1, N))
# state-dependent Gaussian means
mu <- c(1,5)
nobs <- 1000
S <- rep(NA,nobs)
y <- rep(NA,nobs)
# initialise state and observation
S[1] <- sample(1:N, size=1, prob=delta)
y[1] <- rnorm(1, mu[S[1]], 2)
# simulate state and observation processes forward
for(t in 2:nobs) {
  S[t] <- sample(1:N, size=1, prob=Gamma[S[t-1],])
  y[t] <- rnorm(1, mu[S[t]], 2)
}
pal <- c("firebrick","seagreen","navy") # colour palette
plot(y, col=pal[S], type="h")

# Fit and outs 
data.simu<- list(K = 2, y =y,T=nobs,z=S)
stanc("model_SHMM.stan")
fit <- stan(file = 'model_SHMM.stan', data = data.simu,chain=1,cores=3)

traceplot(fit,pars=c("alphas", "betas", "sigmas",'theta'))
outs <- rstan::extract(fit, permuted = TRUE) # return a list of arrays 

#############################
##    From SHMM. AR(1)     ##
#############################
source('simulation.R')
S=rep(states,times=Stimes)
NAS=which(c(1,diff(S))!=0)


# Fit and outs 
data.simu<- list(K = 5, y = t(Obs)[,1],k=c(1,1,1,2,2),T=dim(Obs)[2],z=S,NAS=NAS,u=Stimes,N=Tc)
stanc("model_SHMM.stan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
fit <- stan(file = 'model_SHMM.stan', data = data.simu,chain=3,cores=3)

traceplot(fit,pars=c("alphas", "betas1",'betas2'))
outs <- rstan::extract(fit, permuted = TRUE) # return a list of arrays 

# Sojourn times 
par(mfrow=c(2,3))
for (i in 1:5)
{
  hist(outs$lambda[,i],main=paste('lambda',as.character(i),sep=' '))
  abline(v=lambda[i],col='red')
}


# tmp parameters
par(mfrow=c(5,5))
par(mar=c(2,2,2,2))

for (i in 1:5)
  {
  for (j in 1:5)
  {hist(outs$theta[,i,j],main=paste('p',as.character(i),as.character(j),sep=''))
  abline(v=tpm[i,j],col='red')
  }
}
