data {
  int K; // autoregresive order
  int<lower=0> N; // total lenght of data
  int<lower=0> M; // number of replications 
  int<lower=0> Npred; // number of ts simulated 
  real y[N]; // total data
  int q1[M]; // initial indexes
  int q2[M]; // final indexes
  int ns[M];//final lengths
}
parameters {
  real alpha;
  real beta[K];
  real sigma;
}
model {
  sigma ~ student_t(3, 0, 1);
  // priors for beta? 
  for (m in 1:M)
  {
  for (n in (K+1):ns[m]) {
    real mu = alpha;
    for (k in 1:K)
      mu += beta[k] * y[q1[m]:q2[m]][n-k];
    y[(q1[m]):(q2[m])][n] ~ normal(mu, sigma);
   
  }
}
}
generated quantities {
vector[(N-M)] rss; //residuals 
vector[Npred] ypred;
int ind1;
int ind2;
for (j in 1:K)
{
ypred[j] = y[j];
}
for (n in (K+1):Npred)
{
real mus=alpha;
  for (k in 1:K)
  {
  mus+=beta[k]*ypred[n-k];
  }
  ypred[n]=normal_rng(mus, sigma);
  }
//Residuals
  //ind1=1;
  //ind2=q2[1]-K;
//  for(m in 1:(M-1))
  //{
  //rss[ind1:ind2]; 
  //y[(q1[m]+K):(q2[m])];
  //ind1=ind2+1;
  //ind2=q2[m+1]-((m+1)*K);
  //}
  
  }
