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
######## AR(k) ########
#######################

## 1.Run stan models

N=length(combined.Acx)
K=seq(from=1,to=3,by=1)
n.s[1]=0
M=length(n.s)-1

q1=cumsum(n.s)[1:(length(n.s)-1)]+1
q2=c(cumsum(n.s)[2:length(n.s)])
ns=n.s[2:length(n.s)]
y = combined.Acx
fit.k=list()

# fit AR(k) for k in K
for (j in 1:length(K))
{
ts_data <- list(N = N, y = combined.Acx,K=K[j],M=M,q1=q1,q2=q2,ns=ns,Npred=100)
stanc("model_ark.stan")
fit.k[[j]] <- stan(file = 'model_ark.stan', data = ts_data,chain=3,cores=3)
}



##################################################################
############################## OUTS  #############################
##################################################################


outs.k=list()

###############
## 2. Chains ##
###############

for (p in 1:length(K))
{
traceplot(fit.k[[p]],pars=c("alpha", "beta", "sigma"))
outs.k[[p]] <- rstan::extract(fit.k[[p]], permuted = TRUE) # return a list of arrays 
}

##############################
## 3. Parameter estimation ###
##############################

print(fit.k[[p]], pars=c("alpha", "beta", "sigma", "lp__"), probs=c(.1,.5,.9))
plot(fit.k[[p]],pars=c("alpha", "beta", "sigma"))

#####################
## 3.1. Posteriors ##
#####################

hist(outs.k[[p]]$alpha)
hist(outs.k[[p]]$beta[,1])
hist(outs.k[[p]]$beta[,2])
hist(outs.k$beta[,3])
hist(outs.k[[p]]$sigma)

#####################
## 4. Seasonality? ##
#####################
#---> mirar si los betas quedan dentro de la region Cp (si cumplen polinomio)


###################
## 5. Residuals. ##
###################

## 5.1 Calculation in R 

residuals=list()
  
# usamos la media de las posteriores
for (j in 1:length(K))
{
  residuals[[j]]=numeric(length(ts_data$y)-M*K[j])
  beta.p=apply(outs.k[[j]]$beta,MARGIN=2,FUN=mean)
  alpha.p=mean(outs.k[[j]]$alpha)
  ind=1
  for (m in 1:M)
  {    
  for (n in (K[j]):(ts.lengths[m]-1))
    {
         residuals[[j]][ind]=y[q1[m]+n]-(alpha.p+t(beta.p)%*%y[(q1[m]+n-K[j]):(q1[m]+n-1)])
        ind=ind+1
         }
  
  }
}


#n.posterior=length(outs.k$alpha)
#which=sample(seq(1,n.posterior),size=n.posterior)
#residuals=matrix(NA,ncol=n.posterior,nrow=length(ts_data$y)-M*K)


# for (j in 1:n.posterior)
# {
#   ind=1
#   for(m in 1:M)
#   {
#   for (n in K:(ts.lengths[m]-1))
#   {
#     residuals[ind,j]=y[q1[m]+n]-(alpha.sam[j]+t(beta.sam[j,])%*%y[(q1[m]+n-K):(q1[m]+n-1)])
#     ind=ind+1
#     }
#     
#   }  
# } 

### Plot standarized residual
stand.residuals=list()
for (j in 1:length(K))
{
stand.residuals[[j]]=residuals[[j]]-mean (residuals[[j]])/(sd(residuals[[j]])) 

}

par(mfrow=c(length(K),1))
for (k in 1:length(K))
{
  plot(stand.residuals[[k]],type='l', main= paste('Standardized Residuals AR(', as.character(k),')',sep=''))
}



## 5.2 check independency and normality

par(mfrow=c(length(K),2))
for (j in 1:length(K))
{
  hist(residuals[[j]],main=paste('Residuals AR(',as.character(j),')',sep=''),
       xlab='residuals')
  qqnorm(residuals[[j]])
  qqline(residuals[[j]])
}


# shapiro test for normality
s.test=list()
for (j in 1:length(K))
{
s.test[[j]]=shapiro.test(residuals[[j]]) 
}

s.test
# bajos, rechazamos H0


## 5.3 independency

par(mfrow=c(length(K),3))
for (j in 1:length(K))
{
  ts.plot(as.ts(residuals[[j]]),ylab=paste('Residuals AR(',as.character(j),')',sep=''))
  acf(residuals[[j]],main=paste('Residuals AR(',as.character(j),')',sep=''))
pacf(residuals[[j]],main=paste('Residuals AR(',as.character(j),')',sep=''))
}


# residuals autocorrelation values of AR(1)
signif(acf(residuals[[1]],plot=F)$acf[1:6],2)

## El AR(1) parecería ser el que tiene los residuos menos autocorrelacionados

####################################
## 6. Posterior predictive checks ##
####################################

#ypreds: one for each sample of the posterior distribution 
# Example: 
par(mfrow=c(3,3))



plot(outs.k[[1]]$ypred[i,],type='l')

##############
## 7. WAIC  ##
##############

# 7.1 calculation of lpd (Estimated log pointwise predictive density)
log_lik=list()
lik=list()

# calculation of the log likelihood
rm(i,j,k,n,s,m)
for (k in 1:length(K))
{
  log_lik[[k]]=matrix(NA,nrow=length(ts_data$y)-M*K[k],ncol=length(outs.k[[k]]$alpha))
  lik[[k]]=matrix(NA,nrow=length(ts_data$y)-M*K[k],ncol=length(outs.k[[k]]$alpha))
  
  for (s in 1:length(outs.k[[k]]$alpha))
  {
  beta.s=outs.k[[k]]$beta[s,]
  alpha.s=outs.k[[k]]$alpha[s]
  sigma.s=outs.k[[k]]$sigma[s]
    
  ind=1
  for (m in 1:M) # for each replicate
  {
      
      for (i in K[k]:(ts.lengths[m]-1))
      {
        log_lik[[k]][ind,s]=dnorm(y[q1[m]+i],
                              mean=(alpha.s+t(beta.s)%*%y[(q1[m]+i-K[k]):(q1[m]+i-1)]),
                              sd=sigma.s,log=TRUE)
        
        
      lik[[k]][ind,s]=dnorm(y[q1[m]+i],
                                 mean=(alpha.s+t(beta.s)%*%y[(q1[m]+i-K[k]):(q1[m]+i-1)]),
                                 sd=sigma.s,log=FALSE)
      ind=ind+1  
      }
      
      }
  }
}
  
  
p.waic=rep(NA,length(K))
lpd.waic=rep(NA,length(K))
waic=rep(NA,length(K))
for (k in 1:length(K))
{
lpd.waic[k]=sum(log(apply(lik[[k]],MARGIN=1,FUN=mean)))

log.mean.s=apply(log_lik[[k]],MARGIN=1,FUN=mean)
rm(i,s)
su=rep(0,(length(ts_data$y)-M*K[k]))
for (s in 1:length(outs.k[[k]]$alpha))
{
  for (i in 1:(length(ts_data$y)-M*K[k]))
  {
    su[i]=su[i]+log_lik[[k]][i,s]-log.mean.s[i]  
  }
}

p.waic[k]=sum(su)/(length(outs.k[[k]]$alpha)-1)

waic[k]=-2*(lpd.waic[k]-p.waic[k])
}  

waic


## waic nos da mas majo para AR(1). Nos quedariamos con ese