library(tidyverse)
# This script performes some simultions of the model presented
#######################################################################
############################### Simulations ###########################
#######################################################################

### 1. Set model's parameters

# Sojourn time parameters
lambda=c(20,13,10.8,15.8,17.2) 
M=length(lambda) # Number of possible states
                                                                    
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

betas2=list() 
betas2[[1]]=c(NA,0.1,NA) # beta2 for state 1 for the three acc axis
betas2[[2]]=c(NA,0.2,.35) # beta2 for state 2 for the three acc axis
betas2[[3]]=c(NA,NA,NA) # beta2 for state 3 for the three acc axis
betas2[[4]]=c(0.1,0.1,NA) # beta2 for state 4 for the three acc axis
betas2[[5]]=c(0.15,NA,NA) # beta2 for state 5 for the three acc axis

betas3=list() 
betas3[[1]]=c(NA,NA,NA) # beta3 for state 1 for the three acc axis
betas3[[2]]=c(NA,NA,.05) # beta3 for state 2 for the three acc axis
betas3[[3]]=c(NA,NA,NA) # beta3 for state 3 for the three acc axis
betas3[[4]]=c(NA,0.01,NA) # beta3 for state 4 for the three acc axis
betas3[[5]]=c(NA,NA,NA) # beta3 for state 5 for the three acc axis

k=list(c(1,2,1),c(1,2,3),c(1,1,1),c(2,3,1),c(2,1,1)) # order of each model

# Initial distributions
delta=rep(1,M)/M

# Trasition probability matrix 
tpm <- matrix(c(0.75, 0.14, 0.005, 0.005, 0.1,
                 0.1, 0.7, 0.05, 0.05, 0.1,
                 0.01, 0.15,0.6, 0.09,0.15,
                 0.02, 0.01, 0.01, 0.66,0.3,
                 0.25,0.05,0.20,0.4,0.1), byrow = T, nrow = 5)

# Number of total NAS sequence (i.e, number of changes)
Tc=30

### 2. Simulation of the NAS sequence 
states=numeric(Tc)
for (j in 1:Tc) {
  if (j == 1) {
    states[1] <- sample(x = 1:M, size = 1, prob = delta)
    } else {
      states[j] <- sample(x = 1:M, size = 1, prob = tpm[states[j -1], ])
    }
  }

### 3. Simultation of the sojourn times 
Stimes=numeric(Tc)
for (j in 1:Tc) {
    Stimes[j] <- rpois(1,lambda [states[j]])
}
totL=sum(Stimes)

### 4. Simultation of the observations (AR models)
Obs=matrix(NA,nrow=3,ncol=totL)
count=1 # to index the Obs matrix
for (i in 1:Tc)# change in state 
{
 act.state=states[i]
    for (j in 1:3)# acc axis
   {   
     # set the order of the actual model
     ars=c(alphas[[act.state]][j],betas1[[act.state]][j],
     betas2[[act.state]][j],betas3[[act.state]][j])

     # Simulate the ar mode with arima.sim function
     simu=arima.sim(list(ar=na.omit(ars)), n=Stimes[i])
     Obs[j,(count:(count+length(simu)-1))]=simu
    }
 count=count+length(simu)
}

### 7. Creation of the Data Frame with the simulated data
States=as.factor(rep(rep(states,times=Stimes),3))
Acc=c(Obs[1,],Obs[2,],Obs[3,])
Axis=c(rep('Accx',length(Obs[1,])),rep('Accy',length(Obs[2,])),
       rep('Accz',length(Obs[3,])))
Times=rep(seq(1,totL),3)
Obs.df=data.frame(Times,Acc,Axis,States)

### 6. Plot the simulations
# Plot state sequencies
ggplot(Obs.df %>% filter(Axis=='Accx')) + geom_point(aes(x=Times,y=States,col=States))

# Plot the Observations
ggplot(Obs.df) + geom_line(aes(x=Times,y=Acc,col=Axis),alpha=0.6)+theme_bw()

