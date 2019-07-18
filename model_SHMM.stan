data {
  int<lower=1> K;  // num categories
  int <lower=1> nsims;  // num of time series simulated
  int<lower=0> T[nsims];  // num instances
  int<lower=0> N[nsims];  // num of NAS
  int<lower=0,upper=K> z[nsims,max(T)]; // behaviours
  int<lower=0> u[nsims,max(N)]; // sojuorn times
  int NAS[nsims,max(N)]; // indexes of NAS points
  real y[nsims,3,max(T)];// obs data
  //real y[T];// obs data
  int <lower=0,upper=3> k[3, 5];//AR() orders
 //int <lower=0,upper=2> kz[K]; //ACCX AR() orders
  }
parameters {
  //sojourn times parameters
  real <lower=0> lambda[K];  
  //autorregresive models. Assume AR(1) for now
  real betas1[3,K];
  real betas2[3,K];
  real betas3[3,K];
  real alphas[3,K];
  real<lower=0> sigmas[3,K];
  // K x K tpm
  simplex[K] theta[K]; 
}

model {
  vector[K] log_theta_tr[K];
  real lp;
  real lp_p1;
  real lp_p2;
  
  // prior for alphas and betas ??
  //betas1 ~ student_t(3, 0, 1);
  //betas2 ~ student_t(3, 0, 1);
  //alphas ~ student_t(3, 0, 1);
  //sigmas ~ cauchy(0,2);
  
  
// transpose the tpm and take natural log of entries
for (n_from in 1:K)
  for (n in 1:K)
    log_theta_tr[n, n_from] = log(theta[n_from, n]);
// Compute CDL

for (n in 1:nsims)
{
int count=1;


//First obs
lp = log(0.5) + poisson_lpmf(u[n,1]|lambda[z[n,1]])+normal_lpdf(y[n,1,1] | 0, sigmas[1,z[n,1]])+
normal_lpdf(y[n,1,2] | 0, sigmas[2,z[n,1]])+normal_lpdf(y[n,1,3] | 0, sigmas[3,z[n,1]]);

for (t in 4:(T[n])) { // looping over all the times pdfs
//Acx
if (k[1,z[n,t]]==1 && k[1,z[n,t]]==k[1,z[n,t-1]])//AR(1)
  {lp_p1 =lp + normal_lpdf(y[n,1,t] |(alphas[1,z[n,t]]+betas1[1,z[n,t]]*y[n,1,t-1]), sigmas[1,z[n,t]]);}
else if (k[1,z[n,t]]==2 && k[1,z[n,t]]==k[1,z[n,t-1]]&& k[1,z[n,t]]==k[1,z[n,t-2]])//AR(2)
{lp_p1 =lp + normal_lpdf(y[n,1,t] |(alphas[1,z[n,t]]+betas1[1,z[n,t]]*y[n,1,t-1]+betas2[1,z[n,t]]*y[n,1,t-2]), sigmas[1,z[n,t]]);}

else if (k[1,z[n,t]]==3 && k[1,z[n,t]]==k[1,z[n,t-1]]&& k[1,z[n,t]]==k[1,z[n,t-2]]&& k[1,z[n,t]]==k[1,z[n,t-3]])//AR(3)
{lp_p1 =lp + normal_lpdf(y[n,1,t] |(alphas[1,z[n,t]]+betas1[1,z[n,t]]*y[n,1,t-1]+betas2[1,z[n,t]]*y[n,1,t-2]+betas3[1,z[n,t]]*y[n,1,t-3]), sigmas[1,z[n,t]]);}
lp = lp_p1;

//Acy
if (k[2,z[n,t]]==1 && k[2,z[n,t]]==k[2,z[n,t-1]])//AR(1)
 {lp_p1 =lp + normal_lpdf(y[n,2,t] |(alphas[2,z[n,t]]+betas1[2,z[n,t]]*y[n,2,t-1]), sigmas[2,z[n,t]]);
  }
else if (k[2,z[n,t]]==2 && k[2,z[n,t]]==k[2,z[n,t-1]]&& k[2,z[n,t]]==k[2,z[n,t-2]])//AR(2)
{lp_p1 =lp + normal_lpdf(y[n,2,t] |(alphas[2,z[n,t]]+betas1[2,z[n,t]]*y[n,2,t-1]+betas2[2,z[n,t]]*y[n,2,t-2]), sigmas[1,z[n,t]]);}

else if (k[2,z[n,t]]==3 && k[2,z[n,t]]==k[2,z[n,t-1]]&& k[2,z[n,t]]==k[2,z[n,t-2]]&& k[2,z[n,t]]==k[2,z[n,t-3]])//AR(3)
{lp_p1 =lp + normal_lpdf(y[n,2,t] |(alphas[2,z[n,t]]+betas1[2,z[n,t]]*y[n,2,t-1]+betas2[2,z[n,t]]*y[n,2,t-2]+betas3[2,z[n,t]]*y[n,2,t-3]), sigmas[2,z[n,t]]);}
lp = lp_p1;


//Acz
if (k[3,z[n,t]]==1 && k[3,z[n,t]]==k[3,z[n,t-1]])//AR(1)
  {lp_p1 =lp + normal_lpdf(y[n,3,t] |(alphas[3,z[n,t]]+betas1[3,z[n,t]]*y[n,3,t-1]), sigmas[3,z[n,t]]);}
  
else if (k[3,z[n,t]]==2 && k[3,z[n,t]]==k[3,z[n,t-1]]&& k[3,z[n,t]]==k[3,z[n,t-2]])//AR(2)
{lp_p1 =lp + normal_lpdf(y[n,3,t] |(alphas[3,z[n,t]]+betas1[3,z[n,t]]*y[n,3,t-1]+betas2[3,z[n,t]]*y[n,3,t-2]), sigmas[3,z[n,t]]);}

else if (k[3,z[n,t]]==3 && k[3,z[n,t]]==k[3,z[n,t-1]]&& k[3,z[n,t]]==k[3,z[n,t-2]]&& k[3,z[n,t]]==k[3,z[n,t-3]])//AR(3)
{lp_p1 =lp + normal_lpdf(y[n,3,t] |(alphas[3,z[n,t]]+betas1[3,z[n,t]]*y[n,3,t-1]+betas2[3,z[n,t]]*y[n,3,t-2]+betas3[3,z[n,t]]*y[n,3,t-3]), sigmas[3,z[n,t]]);}
lp = lp_p1;
}

for (m in NAS[n,2:N[n]]) { // looping over NAS
  lp_p2 =lp +  log_theta_tr[z[n,m-1],z[n,m]] + poisson_lpmf(u[n,count]|lambda[z[n,m]]);
  lp = lp_p2;
  count=count+1;
}
}
target += lp;
}

