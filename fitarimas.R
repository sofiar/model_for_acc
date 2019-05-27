#### Modelo autorregresivo
library(ggplot2)
library(gridExtra)
library(dplyr)
library(marima)
source('patterns.R')


### Beahaviour: Eating


tss=list()
for (i in 1:length(Eating))
{
ts=Eating[[i]] %>% select(Acc_x,Acc_y,Acc_z,VeDBA,Pitch.angle,DateTime)

# Windows summarize every (1/5 min)

w=floor(40/5)
n=length(ts$Acc_x)
group=rep(seq(1,floor(n/w)),each=w)
group=c(group,rep(floor(n/w)+1,n-length(group)))

ts$group=as.factor(group)
tss[[i]]=ts %>% group_by(group) %>% summarise(Acx=mean(Acc_x),Acy=mean(Acc_y),
Acz=mean(Acc_z),Pitch=mean(Pitch.angle),VeDBA=mean(VeDBA))

# plot(ts$Acc_x,type='l')
# plot(ts.s$Acx,col='red',type='l')

}

## Combine Time Series 
# combined.Acx=c()
# combined.Acy=c()
# combined.Acz=c()
# for (i in 1: length(Eating))
# {
#   combined.Acx=c(combined.Acx,c(tss[[i]] %>% select(Acx))[[1]],rep(NA,200))
#   combined.Acy=c(combined.Acy,c(tss[[i]] %>% select(Acy))[[1]],rep(NA,200))
#   combined.Acz=c(combined.Acz,c(tss[[i]] %>% select(Acz))[[1]],rep(NA,200))
# }
# 
# combined.Acx=ts(combined.Acx)
# combined.Acy=ts(combined.Acy)
# combined.Acz=ts(combined.Acz)
# 
# plot(combined.Acx)
# lines(combined.Acy,col='red')
# lines(combined.Acz,col='blue')

library(forecast)
auto.arima(combined.Acx)

# autocorrelation funcion 

ts.s=tss[[3]]
acf(ts.s$Acx)
pacf(ts.s$Acx)
p=10 ## 2 seconds

## Fit with MLE (Yule-Walker ?)
ts.fit=arima(ts.s$Acx, order = c(p,0,0),method="ML")
ts.fit$coef
#mle of the innovations variance i.e the residual sum of squares divided by n-p
ts.fit$sigma2
# aic value
ts.fit$aic
# Segun el pradp
nn=ts.fit$nobs
aicc=nn*(1+p/nn)/(1-(p+2)/nn)+nn*log(ts.fit$sigma2)
bic=log(nn)*p+nn*log(ts.fit$sigma2)

###########################################################################
#####################   Diagnosis- Residuals  #############################
###########################################################################

###############
# Correlation # 
###############

acf(ts.fit$residuals)
pacf(ts.fit$residuals)

# Box-Ljung
library(FitAR)
boxresult=LjungBoxTest (ts.fit$residuals,k=5,StartLag=1)
#pvalues
boxresult[,3]

#############
# Cero Mean #
#############
a.r=mean(ts.fit$residuals)
a.sd=sd(ts.fit$residuals)
val=a.r/(sqrt(a.sd/length(ts.fit$residuals)))
pnorm(val)

#################
# Homocedasticy #
#################
plot(ts.fit$residuals)



#############
# Normality #
#############

qqnorm(ts.fit$residuals)
qqline(ts.fit$residuals)
shapiro.test(ts.fit$residuals) # bajo


#### Marima models. Capitulo 19 Peña
library(marima)
