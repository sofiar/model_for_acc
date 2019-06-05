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
