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
plot(combined.Acy,type='l')

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
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

N=length(combined.Acx)
K=seq(from=1,to=3,by=1)
n.s[1]=0
M=length(n.s)-1

q1=cumsum(n.s)[1:(length(n.s)-1)]+1
q2=c(cumsum(n.s)[2:length(n.s)])
ns=n.s[2:length(n.s)]
yx = combined.Acx
yy = combined.Acy
yz = combined.Acz

fitx.k=list()
fity.k=list()
fitz.k=list()



# fit AR(k) for Acx
for (j in 1:length(K))
{
ts_dataX <- list(N = N, y = combined.Acx,K=K[j],M=M,q1=q1,q2=q2,ns=ns,Npred=100)
stanc("model_ark.stan")
fitx.k[[j]] <- stan(file = 'model_ark.stan', data = ts_dataX,chain=3,cores=3)
}

# fit AR(k) for Acy
for (j in 1:length(K))
{
  ts_dataY <- list(N = N, y = combined.Acy,K=K[j],M=M,q1=q1,q2=q2,ns=ns,Npred=100)
  stanc("model_ark.stan")
  fity.k[[j]] <- stan(file = 'model_ark.stan', data = ts_dataY,chain=3,cores=3)
}


for (j in 1:length(K))
{
  ts_dataZ <- list(N = N, y = combined.Acz,K=K[j],M=M,q1=q1,q2=q2,ns=ns,Npred=100)
  stanc("model_ark.stan")
  fitz.k[[j]] <- stan(file = 'model_ark.stan', data = ts_dataZ,chain=3,cores=3)
}

##################################################################
############################## OUTS  #############################
##################################################################


outsx.k=list()
outsy.k=list()
outsz.k=list()

###############
## 2. Chains ##
###############

# Acx
for (p in 1:length(K))
{
traceplot(fitx.k[[p]],pars=c("alpha", "beta", "sigma"))
outsx.k[[p]] <- rstan::extract(fitx.k[[p]], permuted = TRUE) # return a list of arrays 
}

# Acy
for (p in 1:length(K))
{
  traceplot(fity.k[[p]],pars=c("alpha", "beta", "sigma"))
  outsy.k[[p]] <- rstan::extract(fity.k[[p]], permuted = TRUE) # return a list of arrays 
}

# Acz
for (p in 1:length(K))
{
  traceplot(fitz.k[[p]],pars=c("alpha", "beta", "sigma"))
  outsz.k[[p]] <- rstan::extract(fitz.k[[p]], permuted = TRUE) # return a list of arrays 
}

##############################
## 3. Parameter estimation ###
##############################

print(fitx.k[[p]], pars=c("alpha", "beta", "sigma", "lp__"), probs=c(.1,.5,.9))
plot(fitx.k[[p]],pars=c("alpha", "beta", "sigma"))

print(fity.k[[p]], pars=c("alpha", "beta", "sigma", "lp__"), probs=c(.1,.5,.9))
plot(fity.k[[p]],pars=c("alpha", "beta", "sigma"))

print(fitz.k[[p]], pars=c("alpha", "beta", "sigma", "lp__"), probs=c(.1,.5,.9))
plot(fitz.k[[p]],pars=c("alpha", "beta", "sigma"))

#####################
## 3.1. Posteriors ##
#####################

hist(outsx.k[[p]]$alpha)
hist(outsx.k[[p]]$beta[,1])
hist(outsx.k[[p]]$beta[,2])
hist(outsx.k$beta[,3])
hist(outsx.k[[p]]$sigma)

hist(outsy.k[[p]]$alpha)
hist(outsy.k[[p]]$beta[,1])
hist(outsy.k[[p]]$beta[,2])
hist(outsy.k$beta[,3])
hist(outsy.k[[p]]$sigma)

hist(outsz.k[[p]]$alpha)
hist(outsz.k[[p]]$beta[,1])
hist(outsz.k[[p]]$beta[,2])
hist(outsz.k$beta[,3])
hist(outsz.k[[p]]$sigma)

#####################
## 4. Seasonality? ##
#####################
#---> mirar si los betas quedan dentro de la region Cp (si cumplen polinomio)


###################
## 5. Residuals. ##
###################

## 5.1 Calculation in R 
# usamos la media de las posteriores

residuals.x=list()
residuals.y=list()
residuals.z=list()

# Acx
for (j in 1:length(K))
{
  residuals.x[[j]]=numeric(length(ts_dataX$y)-M*K[j])
  beta.p=apply(outsx.k[[j]]$beta,MARGIN=2,FUN=mean)
  alpha.p=mean(outsx.k[[j]]$alpha)
  ind=1
  for (m in 1:M)
  {    
  for (n in (K[j]):(ts.lengths[m]-1))
    {
         residuals.x[[j]][ind]=yx[q1[m]+n]-(alpha.p+t(beta.p)%*%yx[(q1[m]+n-K[j]):(q1[m]+n-1)])
        ind=ind+1
         }
  
  }
}

# Acy
for (j in 1:length(K))
{
  residuals.y[[j]]=numeric(length(ts_dataY$y)-M*K[j])
  beta.p=apply(outsy.k[[j]]$beta,MARGIN=2,FUN=mean)
  alpha.p=mean(outsy.k[[j]]$alpha)
  ind=1
  for (m in 1:M)
  {    
    for (n in (K[j]):(ts.lengths[m]-1))
    {
      residuals.y[[j]][ind]=yy[q1[m]+n]-(alpha.p+t(beta.p)%*%yy[(q1[m]+n-K[j]):(q1[m]+n-1)])
      ind=ind+1
    }
    
  }
}


# Acz
for (j in 1:length(K))
{
  residuals.z[[j]]=numeric(length(ts_dataZ$y)-M*K[j])
  beta.p=apply(outsz.k[[j]]$beta,MARGIN=2,FUN=mean)
  alpha.p=mean(outsz.k[[j]]$alpha)
  ind=1
  for (m in 1:M)
  {    
    for (n in (K[j]):(ts.lengths[m]-1))
    {
      residuals.z[[j]][ind]=yz[q1[m]+n]-(alpha.p+t(beta.p)%*%yz[(q1[m]+n-K[j]):(q1[m]+n-1)])
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
stand.residualsX=list()
stand.residualsY=list()
stand.residualsZ=list()
for (j in 1:length(K))
{
stand.residualsX[[j]]=residuals.x[[j]]-mean (residuals.x[[j]])/(sd(residuals.x[[j]]))
stand.residualsY[[j]]=residuals.y[[j]]-mean (residuals.y[[j]])/(sd(residuals.y[[j]]))
stand.residualsZ[[j]]=residuals.z[[j]]-mean (residuals.z[[j]])/(sd(residuals.z[[j]]))

}

par(mfrow=c(length(K),1))
for (k in 1:length(K))
{
  plot(stand.residualsZ[[k]],type='l', main= paste('Standardized Residuals AR(', as.character(k),')',sep=''))
}



## 5.2 check independency and normality

# Acx
par(mfrow=c(1,2))
for (j in 1:length(K))
{
  hist(residuals.x[[j]],main=paste('Residuals AR(',as.character(j),')',sep=''),
       xlab='residuals')
  qqnorm(residuals.x[[j]])
  qqline(residuals.x[[j]])
}


# Acy
par(mfrow=c(1,2))
for (j in 1:length(K))
{
  hist(residuals.y[[j]],main=paste('Residuals AR(',as.character(j),')',sep=''),
       xlab='residuals')
  qqnorm(residuals.y[[j]])
  qqline(residuals.y[[j]])
}

#Acz
par(mfrow=c(1,2))
for (j in 1:length(K))
{
  hist(residuals.z[[j]],main=paste('Residuals AR(',as.character(j),')',sep=''),
       xlab='residuals')
  qqnorm(residuals.z[[j]])
  qqline(residuals.z[[j]])
}


# shapiro test for normality
s.test.x=list()
s.test.y=list()
s.test.z=list()
for (j in 1:length(K))
{
s.test.x[[j]]=shapiro.test(residuals.x[[j]]) 
s.test.y[[j]]=shapiro.test(residuals.y[[j]]) 
s.test.z[[j]]=shapiro.test(residuals.z[[j]]) 
}

s.test.x
s.test.y
s.test.z
# bajos, rechazamos H0


## 5.3 independency

#Acx
par(mfrow=c(length(K),3))
for (j in 1:length(K))
{
  ts.plot(as.ts(residuals.x[[j]]),ylab=paste('Residuals AR(',as.character(j),')',sep=''))
  acf(residuals.x[[j]],main=paste('Residuals AR(',as.character(j),')',sep=''))
 pacf(residuals.x[[j]],main=paste('Residuals AR(',as.character(j),')',sep=''))
}

#Acy
par(mfrow=c(length(K),3))
for (j in 1:length(K))
{
  ts.plot(as.ts(residuals.y[[j]]),ylab=paste('Residuals AR(',as.character(j),')',sep=''))
  acf(residuals.y[[j]],main=paste('Residuals AR(',as.character(j),')',sep=''))
  pacf(residuals.y[[j]],main=paste('Residuals AR(',as.character(j),')',sep=''))
}

#Acz
par(mfrow=c(length(K),3))
for (j in 1:length(K))
{
  ts.plot(as.ts(residuals.z[[j]]),ylab=paste('Residuals AR(',as.character(j),')',sep=''))
  acf(residuals.z[[j]],main=paste('Residuals AR(',as.character(j),')',sep=''))
  pacf(residuals.z[[j]],main=paste('Residuals AR(',as.character(j),')',sep=''))
}

# residuals autocorrelation values of AR(1)
signif(acf(residuals[[1]],plot=F)$acf[1:6],2)

## En Acx AR(1) parecería ser el que tiene los residuos menos autocorrelacionados

####################################
## 6. Posterior predictive checks ##
####################################
#ypreds: one for each sample of the posterior distribution 


## 6.1 First visual analysis 

par(mfrow=c(2,1))
for (i in 1:2)
{
plot(outsz.k[[p]]$ypred[i,5:100],type='l')
}

# Real ts: Acx 
par(mfrow=c(2,1))
plot(ts(eating_tss[[5]] %>% select(Acx))[1:100],type='l')
plot(ts(eating_tss[[3]] %>% select(Acx))[1:100],type='l')
plot(ts(eating_tss[[4]] %>% select(Acx)),type='l',xlim=c(0,100))

# Real ts: Acy
par(mfrow=c(2,1))
plot(ts(eating_tss[[5]] %>% select(Acy))[1:100],type='l')
plot(ts(eating_tss[[3]] %>% select(Acy))[1:100],type='l')
plot(ts(eating_tss[[4]] %>% select(Acy)),type='l',xlim=c(0,100))

# Real ts: Real Acz 
par(mfrow=c(2,1))
plot(ts(eating_tss[[5]] %>% select(Acz))[1:100],type='l')
plot(ts(eating_tss[[3]] %>% select(Acz))[1:100],type='l')
plot(ts(eating_tss[[4]] %>% select(Acz)),type='l',xlim=c(0,100))


## 6.2 Compute summary statistics

# WATCH OUT: some summary statistics could need that the simulated traj have the
# same length of the raeal ts


##### TEST1: frequency of “switches” (gelman p. 165)

test1 <- function (y){
  n <- length (y)
  y.lag <- c (NA, y[1:(n-1)])
  y.lag2 <- c (NA, NA, y[1:(n-2)])
  sum (sign(y-y.lag) != sign(y.lag-y.lag2), na.rm=TRUE)
}

# AcX
test.repX <- rep (NA, 3000)
for (s in 1:3000){
  test.repX[s] <- test1(outsx.k[[p]]$ypred[s,5:100])
}

hist(test.repX)
abline(v=test1(ts(eating_tss[[5]] %>% select(Acx))[1:95]),col='red')
abline(v=test1(ts(eating_tss[[3]] %>% select(Acx))[1:95]),col='red')

# AcY
test.repY <- rep (NA, 3000)
for (s in 1:3000){
  test.repY[s] <- test1(outsy.k[[p]]$ypred[s,5:100])
}

hist(test.repY)
abline(v=test1(ts(eating_tss[[5]] %>% select(Acy))[1:95]),col='red')
abline(v=test1(ts(eating_tss[[3]] %>% select(Acy))[1:95]),col='red')

# AcZ
test.repZ <- rep (NA, 3000)
for (s in 1:3000){
  test.repZ[s] <- test1(outsz.k[[p]]$ypred[s,5:100])
}

hist(test.repZ)
abline(v=test1(ts(eating_tss[[5]] %>% select(Acz))[1:95]),col='red')
abline(v=test1(ts(eating_tss[[3]] %>% select(Acz))[1:95]),col='red')

##### TEST2: Mean time between peaks

pks=function (x, thresh ) 
{
  pks <- which(diff(sign(diff(x))) < 0) + 2
  pks=pks[x[pks - 1] - x[pks] > thresh]
  return(pks)
}


# mean time between pks (thresh)
test2=function (p)
{
  np=length(p)
  if(np>1)
  {
  tbp=numeric(np-1)
  for (i in 2:np)
  {
    tbp[i-1]=p[i]-p[i-1]
  }
  exit=mean(tbp)
  }
else
  {exit=NA}
  return(exit)
  }

#Acx
test2X.rep <- rep (NA, 3000)
for (s in 1:3000){
  test2X.rep[s] <- test2(pks(outsx.k[[p]]$ypred[s,5:100],thresh = 0.1))
}

hist(test2X.rep,main='Mean Time between peaks')
abline(v=test2X(pks(ts(eating_tss[[5]] %>% select(Acx))[1:95],thresh=0.1)),col='red')
abline(v=test2X(pks(ts(eating_tss[[3]] %>% select(Acx))[1:95],thresh=0.1)),col='red')


#Acy
test2Y.rep <- rep (NA, 3000)
for (s in 1:3000){
  test2Y.rep[s] <- test2(pks(outsy.k[[p]]$ypred[s,5:100],thresh = 0.1))
}

hist(test2Y.rep,main='Mean Time between peaks')
abline(v=test2(pks(ts(eating_tss[[5]] %>% select(Acy))[1:95],thresh=0.1)),col='red')
abline(v=test2(pks(ts(eating_tss[[3]] %>% select(Acy))[1:95],thresh=0.1)),col='red')


#Acz
test2Z.rep <- rep (NA, 3000)
for (s in 1:3000){
  test2Z.rep[s] <- test2(pks(outsz.k[[p]]$ypred[s,5:100],thresh = 0.1))
}

hist(test2Z.rep,main='Mean Time between peaks')
abline(v=test2(pks(ts(eating_tss[[5]] %>% select(Acz))[1:95],thresh=0.1)),col='red')
abline(v=test2(pks(ts(eating_tss[[3]] %>% select(Acz))[1:95],thresh=0.1)),col='red')


##### TEST3: Number of peaks (thresh)


#Acx
test3X.rep <- rep (NA, 3000)
for (s in 1:3000){
  test3X.rep[s] <- length(pks(outsx.k[[p]]$ypred[s,5:100],thresh = 0.1))
}

hist(test3X.rep,main='Number of peaks')
abline(v=length(pks(ts(eating_tss[[5]] %>% select(Acx))[1:95],thresh = 0.1)),col='red')
abline(v=length(pks(ts(eating_tss[[3]] %>% select(Acx))[1:95],thresh = 0.1)),col='red')
         
#Acy
test3Y.rep <- rep (NA, 3000)
for (s in 1:3000){
  test3Y.rep[s] <- length(pks(outsy.k[[p]]$ypred[s,5:100],thresh = 0.1))
}

hist(test3Y.rep,main='Number of peaks')
abline(v=length(pks(ts(eating_tss[[5]] %>% select(Acy))[1:95],thresh = 0.1)),col='red')
abline(v=length(pks(ts(eating_tss[[3]] %>% select(Acy))[1:95],thresh = 0.1)),col='red')

#Acz
test3Z.rep <- rep (NA, 3000)
for (s in 1:3000){
  test3Z.rep[s] <- length(pks(outsz.k[[p]]$ypred[s,5:100],thresh = 0.1))
}

hist(test3Z.rep,main='Number of peaks')
abline(v=length(pks(ts(eating_tss[[5]] %>% select(Acz))[1:95],thresh = 0.1)),col='red')
abline(v=length(pks(ts(eating_tss[[3]] %>% select(Acz))[1:95],thresh = 0.1)),col='red')

##### TEST4: Max (and Range Max-Min)

#Acx
test4x.rep <- rep (NA, 3000)
for (s in 1:3000){
  test4x.rep[s] <- max(outsx.k[[p]]$ypred[s,5:100])#-min(outs.k[[1]]$ypred[s,5:100])
}

hist(test4x.rep,main='Maximum value')
abline(v=max(ts(eating_tss[[5]] %>% select(Acx))[1:95]),col='red')
abline(v=max(ts(eating_tss[[3]] %>% select(Acx))[1:95]),col='red')

#Acy
test4y.rep <- rep (NA, 3000)
for (s in 1:3000){
  test4y.rep[s] <- max(outsy.k[[p]]$ypred[s,5:100])#-min(outs.k[[1]]$ypred[s,5:100])
}

hist(test4y.rep,main='Maximum value')
abline(v=max(ts(eating_tss[[5]] %>% select(Acy))[1:95]),col='red')
abline(v=max(ts(eating_tss[[3]] %>% select(Acy))[1:95]),col='red')

#Acz
test4z.rep <- rep (NA, 3000)
for (s in 1:3000){
  test4z.rep[s] <- max(outsz.k[[p]]$ypred[s,5:100])#-min(outs.k[[1]]$ypred[s,5:100])
}

hist(test4z.rep,main='Maximum value')
abline(v=max(ts(eating_tss[[5]] %>% select(Acz))[1:95]),col='red')
abline(v=max(ts(eating_tss[[3]] %>% select(Acz))[1:95]),col='red')

##### TEST5: Mean difference between consecutive values

#Acx
test5x.rep <- rep (NA, 3000)
for (s in 1:3000){
  test5x.rep[s] <- mean(diff(outsx.k[[p]]$ypred[s,5:100]))
}

hist(test5x.rep,main='Mean difference between consecutive values')
abline(v=mean(diff(ts(eating_tss[[5]] %>% select(Acx))[1:95])),col='red')
abline(v=mean(diff(ts(eating_tss[[3]] %>% select(Acx))[1:95])),col='red')

#Acy
test5y.rep <- rep (NA, 3000)
for (s in 1:3000){
  test5y.rep[s] <- mean(diff(outsy.k[[p]]$ypred[s,5:100]))
}

hist(test5y.rep,main='Mean difference between consecutive values')
abline(v=mean(diff(ts(eating_tss[[5]] %>% select(Acy))[1:95])),col='red')
abline(v=mean(diff(ts(eating_tss[[3]] %>% select(Acy))[1:95])),col='red')


#Acy
test5z.rep <- rep (NA, 3000)
for (s in 1:3000){
  test5z.rep[s] <- mean(diff(outsz.k[[p]]$ypred[s,5:100]))
}

hist(test5z.rep,main='Mean difference between consecutive values')
abline(v=mean(diff(ts(eating_tss[[5]] %>% select(Acz))[1:95])),col='red')
abline(v=mean(diff(ts(eating_tss[[3]] %>% select(Acz))[1:95])),col='red')


##### TEST6: Mean standard deviation over a window 


test6=function(x,w)
{
  nn=floor(length(x)/w)
  ssd=numeric(nn)
  for (i in 1:nn)
  {
    ssd[i]=sd(x[(w*(i-1)+1):(w*(i-1)+w)])
  }
return(mean(ssd))
}

w=5

#Acx
test6x.rep <- rep (NA, 3000)
for (s in 1:3000){
  test6x.rep[s] <- test6(outsx.k[[p]]$ypred[s,5:100],w)
}

hist(test6x.rep,main=paste('Mean sd over a window(',as.character(w),')',sep=''))
abline(v=test6(ts(eating_tss[[5]] %>% select(Acx))[1:95],w),col='red')
abline(v=test6(ts(eating_tss[[3]] %>% select(Acx))[1:95],w),col='red')

#Acy
test6y.rep <- rep (NA, 3000)
for (s in 1:3000){
  test6y.rep[s] <- test6(outsy.k[[p]]$ypred[s,5:100],w)
}

hist(test6y.rep,main=paste('Mean sd over a window(',as.character(w),')',sep=''))
abline(v=test6(ts(eating_tss[[5]] %>% select(Acy))[1:95],w),col='red')
abline(v=test6(ts(eating_tss[[3]] %>% select(Acy))[1:95],w),col='red')

#Acz
test6z.rep <- rep (NA, 3000)
for (s in 1:3000){
  test6z.rep[s] <- test6(outsz.k[[p]]$ypred[s,5:100],w)
}

hist(test6z.rep,main=paste('Mean sd over a window(',as.character(w),')',sep=''))
abline(v=test6(ts(eating_tss[[5]] %>% select(Acz))[1:95],w),col='red')
abline(v=test6(ts(eating_tss[[3]] %>% select(Acz))[1:95],w),col='red')



##############
## 7. WAIC  ##
##############

# Acx
# 7.1 calculation of lpd (Estimated log pointwise predictive density)
log_lik.x=list()
lik.x=list()

# calculation of the log likelihood
rm(i,j,k,n,s,m)
for (k in 1:length(K))
{
  log_lik.x[[k]]=matrix(NA,nrow=length(ts_dataX$y)-M*K[k],ncol=length(outsx.k[[k]]$alpha))
  lik.x[[k]]=matrix(NA,nrow=length(ts_dataX$y)-M*K[k],ncol=length(outsx.k[[k]]$alpha))
  
  for (s in 1:length(outsx.k[[k]]$alpha))
  {
  beta.s=outsx.k[[k]]$beta[s,]
  alpha.s=outsx.k[[k]]$alpha[s]
  sigma.s=outsx.k[[k]]$sigma[s]
    
  ind=1
  for (m in 1:M) # for each replicate
  {
      
      for (i in K[k]:(ts.lengths[m]-1))
      {
        log_lik.x[[k]][ind,s]=dnorm(yx[q1[m]+i],
                              mean=(alpha.s+t(beta.s)%*%yx[(q1[m]+i-K[k]):(q1[m]+i-1)]),
                              sd=sigma.s,log=TRUE)
        
        
      lik.x[[k]][ind,s]=dnorm(yx[q1[m]+i],
                                 mean=(alpha.s+t(beta.s)%*%yx[(q1[m]+i-K[k]):(q1[m]+i-1)]),
                                 sd=sigma.s,log=FALSE)
      ind=ind+1  
      }
      
      }
  }
}
  
  
p.waic.x=rep(NA,length(K))
lpd.waic.x=rep(NA,length(K))
waic.x=rep(NA,length(K))
for (k in 1:length(K))
{
lpd.waic.x[k]=sum(log(apply(lik.x[[k]],MARGIN=1,FUN=mean)))

log.mean.s=apply(log_lik.x[[k]],MARGIN=1,FUN=mean)
rm(i,s)
su=rep(0,(length(ts_dataX$y)-M*K[k]))
for (s in 1:length(outsx.k[[k]]$alpha))
{
  for (i in 1:(length(ts_dataX$y)-M*K[k]))
  {
    su[i]=su[i]+log_lik.x[[k]][i,s]-log.mean.s[i]  
  }
}

p.waic.x[k]=sum(su)/(length(outsx.k[[k]]$alpha)-1)

waic.x[k]=-2*(lpd.waic.x[k]-p.waic.x[k])
}  

waic.x

## waic nos da mas bajo para AR(1). Nos quedariamos con ese

# Acy
# 7.1 calculation of lpd (Estimated log pointwise predictive density)
log_lik.y=list()
lik.y=list()

# calculation of the log likelihood
rm(i,j,k,n,s,m)
for (k in 1:length(K))
{
  log_lik.y[[k]]=matrix(NA,nrow=length(ts_dataY$y)-M*K[k],ncol=length(outsy.k[[k]]$alpha))
  lik.y[[k]]=matrix(NA,nrow=length(ts_dataY$y)-M*K[k],ncol=length(outsy.k[[k]]$alpha))
  
  for (s in 1:length(outsy.k[[k]]$alpha))
  {
    beta.s=outsy.k[[k]]$beta[s,]
    alpha.s=outsy.k[[k]]$alpha[s]
    sigma.s=outsy.k[[k]]$sigma[s]
    
    ind=1
    for (m in 1:M) # for each replicate
    {
      
      for (i in K[k]:(ts.lengths[m]-1))
      {
        log_lik.y[[k]][ind,s]=dnorm(yy[q1[m]+i],
                                    mean=(alpha.s+t(beta.s)%*%yy[(q1[m]+i-K[k]):(q1[m]+i-1)]),
                                    sd=sigma.s,log=TRUE)
        
        
        lik.y[[k]][ind,s]=dnorm(yy[q1[m]+i],
                                mean=(alpha.s+t(beta.s)%*%yy[(q1[m]+i-K[k]):(q1[m]+i-1)]),
                                sd=sigma.s,log=FALSE)
        ind=ind+1  
      }
      
    }
  }
}


p.waic.y=rep(NA,length(K))
lpd.waic.y=rep(NA,length(K))
waic.y=rep(NA,length(K))
for (k in 1:length(K))
{
  lpd.waic.y[k]=sum(log(apply(lik.y[[k]],MARGIN=1,FUN=mean)))
  
  log.mean.s=apply(log_lik.y[[k]],MARGIN=1,FUN=mean)
  rm(i,s)
  su=rep(0,(length(ts_dataY$y)-M*K[k]))
  for (s in 1:length(outsy.k[[k]]$alpha))
  {
    for (i in 1:(length(ts_dataY$y)-M*K[k]))
    {
      su[i]=su[i]+log_lik.y[[k]][i,s]-log.mean.s[i]  
    }
  }
  
  p.waic.y[k]=sum(su)/(length(outsy.k[[k]]$alpha)-1)
  
  waic.y[k]=-2*(lpd.waic.y[k]-p.waic.y[k])
}  

waic.y


# Acz
# 7.1 calculation of lpd (Estimated log pointwise predictive density)
log_lik.z=list()
lik.z=list()

# calculation of the log likelihood
rm(i,j,k,n,s,m)
for (k in 1:length(K))
{
  log_lik.z[[k]]=matrix(NA,nrow=length(ts_dataZ$y)-M*K[k],ncol=length(outsz.k[[k]]$alpha))
  lik.z[[k]]=matrix(NA,nrow=length(ts_dataZ$y)-M*K[k],ncol=length(outsz.k[[k]]$alpha))
  
  for (s in 1:length(outsz.k[[k]]$alpha))
  {
    beta.s=outsz.k[[k]]$beta[s,]
    alpha.s=outsz.k[[k]]$alpha[s]
    sigma.s=outsz.k[[k]]$sigma[s]
    
    ind=1
    for (m in 1:M) # for each replicate
    {
      
      for (i in K[k]:(ts.lengths[m]-1))
      {
        log_lik.z[[k]][ind,s]=dnorm(yz[q1[m]+i],
                                    mean=(alpha.s+t(beta.s)%*%yz[(q1[m]+i-K[k]):(q1[m]+i-1)]),
                                    sd=sigma.s,log=TRUE)
        
        
        lik.z[[k]][ind,s]=dnorm(yz[q1[m]+i],
                                mean=(alpha.s+t(beta.s)%*%yz[(q1[m]+i-K[k]):(q1[m]+i-1)]),
                                sd=sigma.s,log=FALSE)
        ind=ind+1  
      }
      
    }
  }
}


p.waic.z=rep(NA,length(K))
lpd.waic.z=rep(NA,length(K))
waic.z=rep(NA,length(K))
for (k in 1:length(K))
{
  lpd.waic.z[k]=sum(log(apply(lik.z[[k]],MARGIN=1,FUN=mean)))
  
  log.mean.s=apply(log_lik.z[[k]],MARGIN=1,FUN=mean)
  rm(i,s)
  su=rep(0,(length(ts_dataZ$y)-M*K[k]))
  for (s in 1:length(outsz.k[[k]]$alpha))
  {
    for (i in 1:(length(ts_dataZ$y)-M*K[k]))
    {
      su[i]=su[i]+log_lik.z[[k]][i,s]-log.mean.s[i]  
    }
  }
  
  p.waic.z[k]=sum(su)/(length(outsz.k[[k]]$alpha)-1)
  
  waic.z[k]=-2*(lpd.waic.z[k]-p.waic.z[k])
}  

waic.z
