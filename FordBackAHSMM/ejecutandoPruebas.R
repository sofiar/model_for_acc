library(Rcpp)


FF=as.vector(F)
NN=N
NNnorm=(Norm)
SStateIn=as.vector(StateIn)
PDF=as.vector(t(pdf))
Dd=as.vector(t(d))
DD=as.vector(t(D))



Rcpp::sourceCpp('Forwrd.cpp')
a=Forwrd(FF,NN,NNnorm,SStateIn,J,tau,M,delta,PDF,Dd,DD,tpm)
Ford=a$F
ST=a$StateIn


dim(Ford)= c(tau,J) # ?Es asi
dim(ST)= c(tau, J) # ?Es asi
ST=t(ST)
Ford=t(Ford)

Ford[,tau]
FB.Result$Forwrd[,tau]

# for (i in 1:length(Ford[1,]))
# {
# print(Ford[1,i])
# print(FB.Result$Forwrd[1,i])
# }


## Pruebo Backward alg

LL = rep(NA, J*tau)
GG = rep(NA, J*tau)
LL1 = rep(NA, J*tau)
FF=a$F
NN=a$N
SStatein=a$StateIn

GT=FB.Result$Backwrd
LT=FB.Result$Gamma

Rcpp::sourceCpp('Backwrd.cpp')
b=Backwrd(GG, LL, LL1, FF, NN, SStatein,J,tau, M, PDF, 
        Dd, DD,tpm)
Back=b$G
Lb=b$L
dim(Back)= c(tau,J) # ?Es asi
dim(Lb)= c(tau,J) # ?Es asi
Back=t(Back)
Lb=t(Lb)

GT[,tau-3]
Back[,tau-3]
Lb[,tau-1]
LT[,tau-1]


