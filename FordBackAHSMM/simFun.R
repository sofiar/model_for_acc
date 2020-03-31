sim.ts.NP<- function(delta, lambdas, tpm, params, n, seed,M){
  ##############################################################  
  
  #### Description:
  #### simulation of a time series with hsmm structure
  #### with poisson sojourn times and normal observations (1-d)  
  
  #### Inputs: 
  #### delta: vector initial probs (dim J)
  #### lambdas: vector parameters of the poisson (J)
  #### tpm: transition prob matrix (JxJ)
  #### params: list (mu=(),sd=())  
  #### seed: double
  #### n: number of AS  
  #### M: double. maximum runlength
  ##############################################################  
  
  source('funAuxhsmm.R')
  set.seed(seed)
  
  J=length(delta)  
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
  RL=matrix(get.d(rd='pois',J,M,param=list(lambda=lambdas)),nrow=J,byrow=TRUE)
  
  
  for (j in 1:n) {
    i=states[j]
    Stimes[j] =sample(c(1:length(RL[i, RL[i,]!=0])), prob=RL[i, RL[i,]!=0], 1)  
  }
  totL=sum(Stimes)
  
  States=rep(states,times=Stimes)
  
  ## 3. Observations (Nomales independientes)
  
  obs <- rep(NA,totL)
  
  # set first obs
  for (j in 1:totL) {
    mu=params$mu[States[j]]
    sig=params$sd[States[j]]
    obs[j] = rnorm(n=1, mean=mu, sd=sig)
  }
  
  out=data.frame("Obs"=obs,"States"=States)  
  
  return(out)
}