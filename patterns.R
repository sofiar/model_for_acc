library(readxl)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)
source("plotdata.R") # load plot functions
orden=c('Acc_x','Acc_y','Acc_z')



# Import Data set 
hv1 <- read_csv("hv1.csv")
hv2 <- read_csv("hv2.csv")
hv3 <- read_csv("hv3.csv")
hv4 <- read_csv("hv4.csv")
hv5 <- read_csv("hv5.csv")

# transform Behaviours as factors
hv1$Behaviours=as.factor(hv1$Behaviours)
hv2$Behaviours=as.factor(hv2$Behaviours)
hv3$Behaviours=as.factor(hv3$Behaviours)
hv4$Behaviours=as.factor(hv4$Behaviours)
hv5$Behaviours=as.factor(hv5$Behaviours)

# Possible Bahaviours: 
# Walk, Vigilance, Vigilance.down, Search, Ruminating, Shake, Eating, Down.Head


# Lets Plot the first hour of the complete time series 

Fast.walk=list(hv1 %>% filter(Behaviours=='Fast.Walk'),
               hv5 %>% filter(Behaviours=='Walk.Fast'))
Searching=list(hv2 %>% filter(Behaviours=='Search'),
               hv3 %>% filter(Behaviours=='Search'),
               hv5 %>% filter(Behaviours=='Search'))

Eating=list(hv1 %>% filter(Behaviours=='Eating'),
            hv2 %>% filter(Behaviours=='Eating'),
            hv3 %>% filter(Behaviours=='Eating'),
            hv4 %>% filter(Behaviours=='Eating'),
            hv5 %>% filter(Behaviours=='Eating'))

Vigilance=list(hv1 %>% filter(Behaviours=='Vigilance'),
            hv2 %>% filter(Behaviours=='Vigilance'),
            hv3 %>% filter(Behaviours=='Vigilance'),
            hv4 %>% filter(Behaviours=='Vigilance'),
            hv5 %>% filter(Behaviours=='Vigilance'))

Vigilance.down=list(hv2 %>% filter(Behaviours=='Vigilance.Down'),
               hv3 %>% filter(Behaviours=='Vigilance.Down'),
               hv4 %>% filter(Behaviours=='Vigilance.Down'),
               hv5 %>% filter(Behaviours=='Vigilance.Down'))

Walk=list(hv1 %>% filter(Behaviours=='Walk'),
          hv2 %>% filter(Behaviours=='Walk'),
          hv3 %>% filter(Behaviours=='Walk'),
          hv5 %>% filter(Behaviours=='Walk'))


Ruminating=list(hv2 %>% filter(Behaviours=='Ruminating'),
                hv3 %>% filter(Behaviours=='Ruminating'),
                hv4 %>% filter(Behaviours=='Ruminating'),
                hv5 %>% filter(Behaviours=='Ruminating'))

Down.head=list(hv2 %>% filter(Behaviours=='Down.Head'))

Shake=list(hv2 %>% filter(Behaviours=='Shake'),
           hv3 %>% filter(Behaviours=='Shake'))

######################################################################################
################################ Plot patterns #######################################
###################################################################################### 

#################
### Vigilance ###
#################
# ctrl+shit+c
# for (i in 1:length(Vigilance))
# {
#   lstt=Vigilance[[i]]
#   h.list=which(diff(lstt$Event.no.)>1)
#   if (length(h.list)>1)
#   {
#     h.list=c(1,h.list)
#   for (j in 1:(length(h.list)-1))
#   {
#     pp=Plot.accData(lstt[h.list[j]:h.list[j+1],],orden) 
#     ggsave(filename=paste(paste('Vigilance',i,j,sep='_'),'.png',sep=''),plot=pp,
#            path='./patterns_files',width=15,height = 8)
#   
#     }
#      
#   }
#   else
#   {
#     1
#   }
#   
# }
# 
# 
# #################
# ### Fast.Walk ###
# #################
# 
# for (i in 1:length(Fast.walk))
# {
#   lstt=Fast.walk[[i]]
#   h.list=which(diff(lstt$Event.no.)>1)
#   if (length(h.list)>1)
#   {
#     h.list=c(1,h.list)
#     for (j in 1:(length(h.list)-1))
#     {
#       pp=Plot.accData(lstt[h.list[j]:h.list[j+1],],orden) 
#       ggsave(filename=paste(paste('Fast_Walk',i,j,sep='_'),'.png',sep=''),plot=pp,
#              path='./patterns_files',width=15,height = 8)
#       
#     }
#     
#   }
#   else
#   {
#     pp=Plot.accData(lstt,orden) 
#     ggsave(filename=paste(paste('Fast_Walk',i,sep='_'),'.png',sep=''),plot=pp,
#            path='./patterns_files',width=15,height = 8)
#   }
#   
# }
# 
# 
# 
# 
# #################
# ####   Walk  ####
# #################
# 
# for (i in 1:length(Walk))
# {
#   lstt=Walk[[i]]
#   h.list=which(diff(lstt$Event.no.)>1)
#   if (length(h.list)>1)
#   {
#     h.list=c(1,h.list)
#     for (j in 1:(length(h.list)-1))
#     {
#       pp=Plot.accData(lstt[h.list[j]:h.list[j+1],],orden) 
#       ggsave(filename=paste(paste('Walk',i,j,sep='_'),'.png',sep=''),plot=pp,
#              path='./patterns_files',width=15,height = 8)
#       
#     }
#     
#   }
#   else
#   {
#     pp=Plot.accData(lstt,orden) 
#     ggsave(filename=paste(paste('Walk',i,sep='_'),'.png',sep=''),plot=pp,
#            path='./patterns_files',width=15,height = 8)
#   }
#   
# }
# 
# 
# #################
# ####Searching####
# #################
# 
# for (i in 1:length(Searching))
# {
#   lstt=Searching[[i]]
#   h.list=which(diff(lstt$Event.no.)>1)
#   if (length(h.list)>1)
#   {
#     h.list=c(1,h.list)
#     for (j in 1:(length(h.list)-1))
#     {
#       pp=Plot.accData(lstt[h.list[j]:h.list[j+1],],orden) 
#       ggsave(filename=paste(paste('Searching',i,j,sep='_'),'.png',sep=''),plot=pp,
#              path='./patterns_files',width=15,height = 8)
#       
#     }
#     
#   }
#   else
#   {
#     pp=Plot.accData(lstt,orden) 
#     ggsave(filename=paste(paste('Searchingd',i,sep='_'),'.png',sep=''),plot=pp,
#            path='./patterns_files',width=15,height = 8)
#   }
#   
# }
# 
