library(tidyverse)
library(rstan)
library(shinystan)

# Behaviour: Vigilance
beha='vigilance'
source('group_ts.R')

################################################################
###################### Unify Time Series #######################
################################################################



n.s=c()
for (i in 1:length(vigilance_tss))
{
  ts.s=vigilance_tss[[i]]
  n.s[i]=length(ts.s$Acx)
}
ts.lengths=n.s
n.s=c(1,n.s)

combined.Acx=c()
combined.Acy=c()
combined.Acz=c()
for (i in 1: length(Vigilance))
{
  combined.Acx=c(combined.Acx,c(vigilance_tss[[i]] %>% select(Acx))[[1]])
  combined.Acy=c(combined.Acy,c(vigilance_tss[[i]] %>% select(Acy))[[1]])
  combined.Acz=c(combined.Acz,c(vigilance_tss[[i]] %>% select(Acz))[[1]])
}

plot(combined.Acx,type='l')
plot(combined.Acy,type='l')
plot(combined.Acz,type='l')




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
