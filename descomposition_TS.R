############################################################
######## Descomposition TS: Exploratory analysis ###########
############################################################

library(forecast)
library(tidyverse)
library(grid)
library(gridExtra)

################# Creation Times series object ###############


TS=Fast.walk[[1]]$Acc_z
TS=hv1$Acc_x

TS=ts(TS,frequency = 40)
Plot.dataclass(hv1,orden)

####### Descomposition: seasonal, trend and remainder ##########

# forecast library
?mstl
d2.ts=mstl(TS,s.window = 40) 
plot(d2.ts)

# stats library
d.ts=stl(TS,s.window ='periodic') # s.window??
plot(d.ts)

d.ts=decompose(TS)
plot(d.ts)

######################### Fourier Analysis #########################



# Periogram: elimination of the tendency ?? 
ssp = spectrum(hv1 %>% select(Acc_x,Acc_y,Acc_z))
## smooth version. span is the number of spikes in the kernel
spectrum(TS,span=5)
spectrum(TS,span=10)
spectrum(TS$Acc_x,span=30)




?#Filter in the frequency domain
library(spectral)
filter.fft(na.omit(TS$Acc_x))

