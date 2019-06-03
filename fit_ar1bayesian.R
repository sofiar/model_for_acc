#######################
######## AR(1) ########
#######################

## Run stan model
N=length(combined.Acx)
M=length(n.s)-1
resdims=cumsum(n.s-1)
resdims[1]=0

ts_dat <- list(N = N, y = combined.Acx,M=M, ydims=cumsum(n.s),rdims=resdims)
stanc("model_ar1.stan")
fit1 <- stan(file = 'model_ar1.stan', data = ts_dat,chain=2,cores=3)

#######
# Out #
#######

## Chains
traceplot(fit1,pars=c("alpha", "beta", "sigma"))
outs <- rstan::extract(fit1, permuted = TRUE) # return a list of arrays 

## Parameter estimation
print(fit1, pars=c("alpha", "beta", "sigma", "lp__"), probs=c(.1,.5,.9))
plot(fit1,pars=c("alpha", "beta", "sigma"))

## Posteriors
hist(la$alpha)
hist(la$beta) # esta dentro del (-1,1)? ---> estacionalidad
hist(la$sigma)

## Residuals
residuals=as.vector(outs$rss)
# cero mean?
a.r=mean(residuals)
a.sd=sd(residuals)
val=a.r/(sqrt(a.sd/length(residuals)))
pnorm(val)

# Normality?
hist(residuals)
qqnorm(residuals)
qqline(residuals)