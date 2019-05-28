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
ts.lengths=n.s
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
fit1 <- stan(file = 'model_ar1.stan', data = ts_dat,chain=2,cores=3)

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

## 1.Run stan model

N=length(combined.Acx)
K=3
n.s[1]=0
M=length(n.s)-1

q1=cumsum(n.s)[1:(length(n.s)-1)]+1
q2=c(cumsum(n.s)[2:length(n.s)])
ns=n.s[2:length(n.s)]
y = combined.Acx

ts_data <- list(N = N, y = combined.Acx,K=K,M=M,q1=q1,q2=q2,ns=ns,Npred=100)

stanc("model_ark.stan")
fitk <- stan(file = 'model_ark.stan', data = ts_data,chain=3,cores=3)


#######
# Out #
#######

## 2. Chains
traceplot(fitk,pars=c("alpha", "beta", "sigma"))
outs.k <- rstan::extract(fitk, permuted = TRUE) # return a list of arrays 

## 3. Parameter estimation
print(fitk, pars=c("alpha", "beta", "sigma", "lp__"), probs=c(.1,.5,.9))
plot(fitk,pars=c("alpha", "beta", "sigma"))

## 3.1. Posteriors
hist(outs.k$alpha)
hist(outs.k$beta[,1])
hist(outs.k$beta[,2])
hist(outs.k$beta[,3])
hist(outs.k$sigma)

## 4. Seasonality? ---> mirar si los betas quedan dentro de la region Cp (si cumplen polinomio)

## 5. Residuals. 
## 5.1 Calculation in R (BIEN ESTO ?)

n.posterior=length(outs.k$alpha)
which=sample(seq(1,n.posterior),size=n.posterior)

beta.sam=outs.k$beta[which,]  
alpha.sam=outs.k$alpha[which]  

residuals=matrix(NA,ncol=n.posterior,nrow=length(ts_data$y)-M*K)

sum(is.na(residuals))
for (j in 1:n.posterior)
{
  ind=1
  for(m in 1:M)
  {
  for (n in K:(ts.lengths[m]-1))
  {
    residuals[ind,j]=y[q1[m]+n]-(alpha.sam[j]+t(beta.sam[j,])%*%y[(q1[m]+n-K):(q1[m]+n-1)])
    ind=ind+1
    }
    
  }  
} 



## 5.2 check independency and normality

nsim=4
hist(residuals[,nsim])
qqnorm(residuals[,nsim])
qqline(residuals[,nsim])
shapiro.test(residuals[,nsim]) # bajo, rechazamos H0


acf(residuals[,nsim])
pacf(residuals[,nsim])



## 6. Posterior predictive checks: to see if the model fit well to the data and compare models

#ypreds: one for each sample of the posterior distribution 
# Example: 
plot(outs.k$ypred[7,],type='l')


## 7. BIC(Prado) and DIC


