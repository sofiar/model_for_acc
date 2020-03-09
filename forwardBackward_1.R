library(matrixStats)
library(RcppArmadillo)
Rcpp::sourceCpp('rcpp_functions.cpp')
source('simu_decod.R')


error.pred=numeric(nrep)

for (i in 1:nrep)
{
log.forwrd=matrix(NA,ncol=n,nrow = M)
log.backwrd=matrix(NA,ncol=n,nrow = M)

### 1.  Forward probabilities
for(m in 1:M)
{
log.forwrd[m,1]=log_allprobs[[i]][1,m]+log_delta[m]
}

# Recurrencia
for (t in 2:n)
{
  for (m in 1:M)
  {
#    log.forwrd[m,t]=log_allprobs[[i]][t,m]+log_sum_expS(log.forwrd[,t-1],t(log_gamma)[,m])# traspuesta de gamma?
    log.forwrd[m,t]=log_allprobs[[i]][t,m]+logSumExp(c(log.forwrd[,t-1]+log_gamma[,m]))
    
    
      }
}


### 2. Backward probabilities
for(m in 1:M)
{
  log.backwrd[m,n]=0
}

# Recurrencia
for (t in (n-1):1)
{
  for (m in 1:M)
  {
    log.backwrd[m,t]=logSumExp(c(log.backwrd[,t+1]+log_allprobs[[i]][t+1,]+
                                   log_gamma[m,] ))
                                  
  }
}

stateprobs=matrix(NA,nrow=n,ncol=M)
lscale = log_sum_expS(log.forwrd[,n], rep(0,M))
for (t in 1:n)
{
  for (j in 1:M)
  {
    stateprobs[t,j]=exp(log.forwrd[j,t] + log.backwrd[j,t] - lscale)
  }
}

## loss function 0-1
preds=numeric(n)
for (t in 1:n)
{
  preds[t]=which.max(stateprobs[t,])  
}

error.pred[i]= sum(states[,i]==preds)/n

}
hist(error.pred)

