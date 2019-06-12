library(ggplot2)
library(gridExtra)
library(dplyr)
library(marima)
source('patterns.R')

############################################################################
############################ Beahaviour: Eating ############################
############################################################################
if (beha=='eating')
{
eating_tss=list()
for (i in 1:length(Eating))
{
  ts=Eating[[i]] %>% select(Acc_x,Acc_y,Acc_z,VeDBA,Pitch.angle,DateTime)
  
  # Windows summarize every (1/2 min)
  
  w=floor(40/2)
  n=length(ts$Acc_x)
  group=rep(seq(1,floor(n/w)),each=w)
  group=c(group,rep(floor(n/w)+1,n-length(group)))
  
  ts$group=as.factor(group)
  eating_tss[[i]]=ts %>% group_by(group) %>% summarise(Acx=mean(Acc_x),Acy=mean(Acc_y),
                                                Acz=mean(Acc_z),Pitch=mean(Pitch.angle),VeDBA=mean(VeDBA))

}

}
############################################################################
############################ Beahaviour: Vigilance #########################
############################################################################
if (beha=='vigilance')
{
vigilance_tss=list()
for (i in 1:length(Vigilance))
{
  ts=Vigilance[[i]] %>% select(Acc_x,Acc_y,Acc_z,VeDBA,Pitch.angle,DateTime)
  
  # Windows summarize every (1/2 min)
  
  w=floor(40/2)
  n=length(ts$Acc_x)
  group=rep(seq(1,floor(n/w)),each=w)
  group=c(group,rep(floor(n/w)+1,n-length(group)))
  
  ts$group=as.factor(group)
  vigilance_tss[[i]]=ts %>% group_by(group) %>% summarise(Acx=mean(Acc_x),Acy=mean(Acc_y),
                                                       Acz=mean(Acc_z),Pitch=mean(Pitch.angle),VeDBA=mean(VeDBA))
  
}

}

############################################################################
############################ Beahaviour: Searching #########################
############################################################################
if (beha=='searching')
{
  searching_tss=list()
  for (i in 1:length(Searching))
  {
    ts=Searching[[i]] %>% select(Acc_x,Acc_y,Acc_z,VeDBA,Pitch.angle,DateTime)
    
    # Windows summarize every (1/2 min)
    
    w=floor(40/2)
    n=length(ts$Acc_x)
    group=rep(seq(1,floor(n/w)),each=w)
    group=c(group,rep(floor(n/w)+1,n-length(group)))
    
    ts$group=as.factor(group)
    searching_tss[[i]]=ts %>% group_by(group) %>% summarise(Acx=mean(Acc_x),Acy=mean(Acc_y),
                                                            Acz=mean(Acc_z),Pitch=mean(Pitch.angle),VeDBA=mean(VeDBA))
    
  }
  
}