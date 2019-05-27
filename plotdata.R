### R plot functions 
library(animalTrack)
library(tidyverse)
library(gridExtra)
library(ggplot2)
library(plotly)

Plot.accData<-function(OvejaDF,orden)
{
  
  OvejaDF$Event.no.=seq(1,length(OvejaDF$Acc_x))
  
  NLS=c(OvejaDF$Acc_x,OvejaDF$Acc_y,OvejaDF$Acc_z)
  ACC=c(rep('Acc_x',length(OvejaDF$Acc_x)),rep('Acc_y',length(OvejaDF$Acc_y)),rep('Acc_z',length(OvejaDF$Acc_z)))
  event.no.=rep(OvejaDF$Event.no.,3)
  
  AccShepPlot=data.frame('Acc.data'=NLS,'ACC'=ACC,'Event.no.'=event.no.)
  
  
  Plot.Acc=ggplot()+geom_line(data=AccShepPlot,aes(x=Event.no.,y=Acc.data,color=ACC),
                              size=0.2,alpha=0.5)+ ylab('Acc')+
    scale_color_manual(values=c("darkgreen", "red", "blue"), 
                       name="",
                       breaks=c('Acc_x', 'Acc_y','Acc_z'),
                       labels=c(orden[1], orden[2], orden[3]))+
    theme(legend.position = 'top',legend.key.size = unit(1, "cm"))
  
  
  Plot.Mag=ggplot(data=OvejaDF)+geom_line(mapping=aes(x=Event.no.,y=Mag_x),
                                          col="#D55E00",size=0.2,alpha=0.7)+
    geom_line(mapping=aes(x=Event.no.,y=Mag_y),col="#56B4E9",size=0.2,alpha=0.7)+
    geom_line(mapping=aes(x=Event.no.,y=Mag_z),col="#009E73",size=0.2,alpha=0.7)+
    ylab('Mag')+xlab('')
  
  Plot.VedBA=ggplot(data=OvejaDF)+geom_line(mapping=aes(x=Event.no.,y=VeDBA),
                                            col="#CC79A7")+ylab('VeDBA')+xlab('')
  Plot.Pitch=ggplot(data=OvejaDF)+geom_line(mapping=aes(x=Event.no.,y=Pitch.angle),
                                            col="#0072B2")+ylab('Pitch')+xlab('')
  
  
  
  PLot=grid.arrange(rbind(ggplotGrob(Plot.Acc),ggplotGrob(Plot.Mag), 
                          ggplotGrob(Plot.VedBA),ggplotGrob(Plot.Pitch),
                          size = "last"))
  return(PLot)
  
}


Plotly.accData<-function(OvejaDF,orden,means=FALSE)
{
  event.no.=rep(OvejaDF$Event.no.,3)
  datetime=rep(OvejaDF$DateTime,3)
  ACC=c(rep(orden[1],length(OvejaDF$Acc_x)),rep(orden[2],length(OvejaDF$Acc_y)),rep(orden[3],length(OvejaDF$Acc_z)))
  MAG=c(rep('Mag_x',length(OvejaDF$Mag_x)),rep('Mag_y',length(OvejaDF$Mag_y)),rep('Mag_z',length(OvejaDF$Mag_z)))
  
  if(means)
  {
    NLS=c(OvejaDF$MeanAcc_x,OvejaDF$MeanAcc_y,OvejaDF$MeanAcc_z)
    MLS=c(OvejaDF$MeanMag_x,OvejaDF$MeanMag_y,OvejaDF$MeanMag_z)
    
  }
  
  else
  {
    NLS=c(OvejaDF$Acc_x,OvejaDF$Acc_y,OvejaDF$Acc_z)
    MLS=c(OvejaDF$Mag_x,OvejaDF$Mag_y,OvejaDF$Mag_z)
  }
  
  
  AccShepPlot=data.frame('Acc.data'=NLS,'ACC'=ACC,
                         'DateTime'=datetime,'event.no'=event.no.)
  
  MagShepPlot=data.frame('MAG.data'=MLS,'MAG'=MAG,
                         'DateTime'=datetime)
  
  Plot.Acc=plot_ly(data=AccShepPlot,x=~DateTime,y=~Acc.data,type='scatter',color=ACC,mode='lines',alpha=0.5)
  
  Plot.Mag=plot_ly(data=MagShepPlot,x=~DateTime,y=~MAG.data,color=MAG,mode='lines',type='scatter',alpha=0.5)
  Plot.VedBA=plot_ly(data=OvejaDF,x=~DateTime,y=~VeDBA, type = 'scatter',mode='lines')
  
  
  Plot.Pitch=plot_ly(data=OvejaDF, x=~DateTime,y=~Pitch.angle,type = 'scatter',mode='lines')%>%
    layout(yaxis = list(title = "Pitch")) %>%
    layout(yaxis = list(range = c(-80, 80)))
  
  subplot(list(Plot.Acc,Plot.Mag,Plot.VedBA,Plot.Pitch), nrows = 4,
          shareX = TRUE)
  
  
}

Plot.dataclass<-function(OvejaDF,orden)
{
  
  OvejaDF$Event.no.=seq(1,length(OvejaDF$Acc_x))
  beha=ggplot(data=OvejaDF)+geom_point(aes(x=Event.no.,y=Behaviours,
                                           color=Behaviours)) + theme(legend.position="none")+xlab('Time')+ theme_bw() 
  
  B.points=which(OvejaDF$Behaviours[-1] != OvejaDF$Behaviours[-length(OvejaDF$Behaviours)])
  NLS=c(OvejaDF$Acc_x,OvejaDF$Acc_y,OvejaDF$Acc_z)
  ACC=c(rep('Acc_x',length(OvejaDF$Acc_x)),rep('Acc_y',length(OvejaDF$Acc_y)),rep('Acc_z',length(OvejaDF$Acc_z)))
  event.no.=rep(OvejaDF$Event.no.,3)
  
  AccShepPlot=data.frame('Acc.data'=NLS,'ACC'=ACC,'Event.no.'=event.no.)
  
  
  Plot.Acc=ggplot()+geom_line(data=AccShepPlot,aes(x=Event.no.,y=Acc.data,color=ACC),
                              size=0.2,alpha=0.5)+ ylab('Acc')+
    scale_color_manual(values=c("green", "red", "blue"), 
                       name="",
                       breaks=c('Acc_x', 'Acc_y','Acc_z'),
                       labels=c(orden[1], orden[2], orden[3]))+
    theme(legend.position = 'top',legend.key.size = unit(1, "cm"))+
    guides(fill = guide_legend(override.aes = list(shape = 22)))+ 
    geom_vline(xintercept = B.points,col="#999999")+ylab('Acc')+xlab('')+ theme_bw()
  
  Plot.VedBA=ggplot(data=OvejaDF)+geom_line(mapping=aes(x=Event.no.,y=VeDBA),
                                            col="#CC79A7")+
    geom_vline(xintercept = B.points,col="#999999")+ylab('VeDBA')+xlab('')+ theme_bw()
  
  Plot.Pitch=ggplot(data=OvejaDF)+geom_line(mapping=aes(x=Event.no.,y=Pitch.angle),
                                            col="#0072B2")+
    geom_vline(xintercept = B.points,col="#999999")+ylab('Pitch')+xlab('')+ theme_bw()

    g1 <- ggplotGrob(Plot.Acc)
  g2 <- ggplotGrob(Plot.VedBA)
  g3 <- ggplotGrob(Plot.Pitch)
  g4 <- ggplotGrob(beha)
  
  colnames(g1) <- paste0(seq_len(ncol(g1)))
  colnames(g2) <- paste0(seq_len(ncol(g2)))
  colnames(g3) <- paste0(seq_len(ncol(g3)))
  colnames(g4) <- paste0(seq_len(ncol(g4)))
  
  grid.draw(gtable_combine(g1, g2,g3,g4, along=2))
  
}

Convert.Data<-function(OvejaDF)
{
  
  OvejaDF2=OvejaDF
  OvejaDF2$Acc_z=-OvejaDF2$Acc_z
  OvejaDF2$Acc_y=-OvejaDF2$Acc_y
  OvejaDF2$Acc_x=OvejaDF2$Acc_x
  
  OvejaDF2$Acc_z.sm=-OvejaDF2$Acc_z.sm
  OvejaDF2$Acc_y.sm=-OvejaDF2$Acc_y.sm
  OvejaDF2$Acc_x.sm=OvejaDF2$Acc_x.sm
  
  #a=pitch(OvejaDF2$Acc_z.sm,OvejaDF2$Acc_y.sm, OvejaDF2$Acc_x.sm)
  a=pitch(OvejaDF2$Acc_x,OvejaDF2$Acc_y, -OvejaDF2$Acc_z)
  OvejaDF2$Pitch.angle=a*(180/pi)# pasamos a grados
  return(OvejaDF2)  
  
}


vline <- function(x = 0,y00,y11) {
  list(
    type = "line", 
    y0 =y00, 
    y1 = max(AccShepPlot$Acc.data), 
    yref = "paper",
    x0 = x, 
    x1 = x, 
    line = list(dash = 'dot',color='grey')
  )
}

Plotly.dataclass<-function(OvejaDF,orden,means=FALSE)
{
  
  OvejaDF$Event.no.=seq(1,length(OvejaDF$Acc_x))
  beha=ggplot(data=OvejaDF)+geom_point(aes(x=Event.no.,y=Behaviours,
                                           color=Behaviours)) + theme(legend.position="none")+xlab('Time')+ theme_bw() 
  
  B.points=which(OvejaDF$Behaviours[-1] != OvejaDF$Behaviours[-length(OvejaDF$Behaviours)])

  event.no.=rep(OvejaDF$Event.no.,3)
  datetime=rep(OvejaDF$DateTime,3)
  ACC=c(rep(orden[1],length(OvejaDF$Acc_x)),rep(orden[2],length(OvejaDF$Acc_y)),rep(orden[3],length(OvejaDF$Acc_z)))
  MAG=c(rep('Mag_x',length(OvejaDF$Mag_x)),rep('Mag_y',length(OvejaDF$Mag_y)),rep('Mag_z',length(OvejaDF$Mag_z)))
  
  if(means)
  {
    NLS=c(OvejaDF$MeanAcc_x,OvejaDF$MeanAcc_y,OvejaDF$MeanAcc_z)
    MLS=c(OvejaDF$MeanMag_x,OvejaDF$MeanMag_y,OvejaDF$MeanMag_z)
    
  }
  
  else
  {
    NLS=c(OvejaDF$Acc_x,OvejaDF$Acc_y,OvejaDF$Acc_z)
    MLS=c(OvejaDF$Mag_x,OvejaDF$Mag_y,OvejaDF$Mag_z)
  }
  
  
  AccShepPlot=data.frame('Acc.data'=NLS,'ACC'=ACC,
                         'DateTime'=datetime,'event.no'=event.no.)
  
  MagShepPlot=data.frame('MAG.data'=MLS,'MAG'=MAG,
                         'DateTime'=datetime)
  
  cortes=list()
  for (i in 1:length(B.points))
  {
    cortes[[i]]=vline(AccShepPlot$DateTime[B.points[i]],min(AccShepPlot$Acc.data),
                      max(AccShepPlot$Acc.data))
  }
  Plot.Acc=plot_ly(data=AccShepPlot,x=~DateTime,y=~Acc.data,type='scatter',color=ACC,mode='lines',alpha=0.5) %>%
    layout(shapes =cortes)
  
  Plot.Mag=plot_ly(data=MagShepPlot,x=~DateTime,y=~MAG.data,color=MAG,mode='lines',type='scatter',alpha=0.5) %>%
    layout(shapes =cortes)
  Plot.VedBA=plot_ly(data=OvejaDF,x=~DateTime,y=~VeDBA, type = 'scatter',mode='lines') %>%
    layout(shapes =cortes)
  
  
  Plot.Pitch=plot_ly(data=OvejaDF, x=~DateTime,y=~Pitch.angle,type = 'scatter',mode='lines')%>%
    layout(yaxis = list(title = "Pitch")) %>%
    layout(yaxis = list(range = c(-80, 80))) %>%
    layout(shapes =cortes)
  
  subplot(list(Plot.Acc,Plot.Mag,Plot.VedBA,Plot.Pitch), nrows = 4,
          shareX = TRUE)
  
  
  
  
  }
