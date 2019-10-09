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
##    From SHMM.           ##
#############################
source('simulation.R')
S=rep(states,times=Stimes)
NAS=which(c(1,diff(S))!=0)
u=numeric(dim(Obs)[2])
u[NAS]=Stimes
# Fit and outs 
data.simu<- list(K = 5, y = Obs,T=dim(Obs)[2],z=S,NAS=NAS,u=u,N=Tc)
stanc("model_SHMM.stan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
fit <- stan(file = 'model_SHMM.stan', data = data.simu,chain=2,cores=3,control = list(max_treedepth = 15))

########################### Results ################################

outs <- rstan::extract(fit, permuted = TRUE) # return a list of arrays 

# Save workspace
save.image("last_fit.RData")
traceplot(fit,pars=c("alphas", "betas1",'betas2'))
fit_summary <- summary(fit,  probs = c(0.025, 0.05, 0.5, 0.95, 0.975))$summary
#op <- par(mfrow = c(1,2))
#hist(fit_summary[,1], main = "R-hat")#
#hist(fit_summary[,9], main = "n-eff" )
#par(op)


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
  {hist(outs$ta[,i,j],main=paste('p',as.character(i),as.character(j),sep=''))
   abline(v=tpm[i,j],col='red')
  }
}

# Autorregressive parameteres
par(mfrow=c(5,3))
par(mar=c(2,2,2,2))
Ac=c('Accx','Accy','Accz')
# Alphas
for (j in 1:5)
{
  for (i in 1:3)
{
      hist(outs$alpha[,i,j],main=paste('alpha',Ac[i],as.character(j),sep=' '))
  abline(v=alphas[[j]][i],col='red')
  }
 }

#Beta1
for (j in 1:5)
{
  for (i in 1:3)
  {
    hist(outs$betas1[,i,j],main=paste('beta1',Ac[i],as.character(j),sep=' '))
    abline(v=betas1[[j]][i],col='red')
    }
}


#Beta2
for (j in 1:5)
{
  for (i in 1:3)
  {
    hist(outs$betas1[,i,j],main=paste('beta2',Ac[i],as.character(j),sep=' '))
    abline(v=betas1[[j]][i],col='red')
  }
}

