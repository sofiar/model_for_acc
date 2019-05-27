## AR(1) with STAN
library(tidyverse)
library(rstan)
library(shinystan)
source('group_ts.R')

# Behaviour: Eating 
################################################################
###################### Unify Time Series #######################
################################################################

n.s=c()
for (i in 1:length(eating_tss))
{
  ts.s=eating_tss[[i]]
  n.s[i]=length(ts.s$Acx)
}
n.s=c(1,n.s)

combined.Acx=c()
combined.Acy=c()
combined.Acz=c()
for (i in 1: length(Eating))
{
  combined.Acx=c(combined.Acx,c(eating_tss[[i]] %>% select(Acx))[[1]])
  combined.Acy=c(combined.Acy,c(eating_tss[[i]] %>% select(Acy))[[1]])
  combined.Acz=c(combined.Acz,c(eating_tss[[i]] %>% select(Acz))[[1]])
}

plot(combined.Acx,type='l')

################################################################
###################### Descriptive Analysis ####################
################################################################
# Autocorrelation function
for (j in 1:length(Eating))
{ acf(eating_tss[[i]])}
 

##################################################################
###################### Inference with Stan  ######################
##################################################################

#######################
######## AR(1) ########
#######################

## Run stan model
N=length(combined.Acx)
M=length(n.s)-1
resdims=cumsum(n.s-1)
resdims[1]=0

ts_dat <- list(N = N, y = combined.Acx,M=M, ydims=cumsum(n.s),rdims=resdims)
stanc("model_ar1.stan")
fit1 <- stan(file = 'model_ar1.stan', data = ts_dat,chain=4,cores=3)

#######
# Out #
#######

## Chains
traceplot(fit1,pars=c("alpha", "beta", "sigma"))
outs <- rstan::extract(fit1, permuted = TRUE) # return a list of arrays 

## Parameter estimation
print(fit1, pars=c("alpha", "beta", "sigma", "lp__"), probs=c(.1,.5,.9))
plot(fit1,pars=c("alpha", "beta", "sigma"))

## Posteriors
hist(la$alpha)
hist(la$beta) # esta dentro del (-1,1)? ---> estacionalidad
hist(la$sigma)

## Residuals
residuals=as.vector(outs$rss)
# cero mean?
a.r=mean(residuals)
a.sd=sd(residuals)
val=a.r/(sqrt(a.sd/length(residuals)))
pnorm(val)

# Normality?
hist(residuals)
qqnorm(residuals)
qqline(residuals)

#######################
######## AR(k) ########
#######################

N=length(combined.Acx)
K=3
n.s[1]=0
M=length(n.s)-1

q1=cumsum(n.s)[1:(length(n.s)-1)]+1
q2=c(cumsum(n.s)[2:length(n.s)])
ns=n.s[2:length(n.s)]
ts_data <- list(N = N, y = combined.Acx,K=K,M=M,q1=q1,q2=q2,ns=ns,Npred=100)

stanc("model_ark.stan")
fitk <- stan(file = 'model_ark.stan', data = ts_data,chain=3,cores=3)


#######
# Out #
#######

## Chains
traceplot(fitk,pars=c("alpha", "beta", "sigma"))
outs.k <- rstan::extract(fitk, permuted = TRUE) # return a list of arrays 

## Parameter estimation
print(fitk, pars=c("alpha", "beta", "sigma", "lp__"), probs=c(.1,.5,.9))
plot(fitk,pars=c("alpha", "beta", "sigma"))

## Posteriors
hist(la$alpha)
hist(la$beta)
hist(la$sigma)

# Estacionalidad ? ---> mirar si los betas quedan dentro de la region Cp (si cumplen polinomio)

# Residuals. Calculation in R (BIEN ESTO ?)

residuals=numeric(length(ts_data$y))

beta.sam=sample(1,la$beta)  
alpha.sam=sample(1,la$alpha)  


ind1=1;
ind2=q2[1]-K;
for (i in 1:M)
{
residuals[ind1:ind2]=y[(q1[m]+K):(q2[m])]-(alpha.sam+beta.sam * y[q1[m]:q2[m]])
}






### ypreds: one for each sample of the posterior distribution 
# Example: 
plot(outs.k$ypred[7,],type='l')

