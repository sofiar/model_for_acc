data {
  int<lower=1> K;  // num categories
  int<lower=0> T;  // num instances
  int<lower=0> N;  // num of NAS
  int<lower=1,upper=K> z[T]; // behaviours
  int<lower=0> u[N]; // sojuorn times
  int NAS[N]; // indexes of NAS points
  matrix[3,T] y;// obs data
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
  int count=1;
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

//First obs
lp = log(0.5) + poisson_lpmf(u[1]|lambda[z[1]])+normal_lpdf(y[1,1] | 0, sigmas[1,z[1]])+
normal_lpdf(y[1,3] | 0, sigmas[2,z[1]])+normal_lpdf(y[1,3] | 0, sigmas[3,z[1]]);

for (t in 3:T) { // looping over all observations pdfs
//Acx
if (k[1,z[t]]==1 && k[1,z[t]]==k[1,z[t-1]])//AR(1)
  {lp_p1 =lp + normal_lpdf(y[1,t] |(alphas[1,z[t]]+betas1[1,z[t]]*y[1,t-1]), sigmas[1,z[t]]);
  }
else if (k[1,z[t]]==2 && k[1,z[t]]==k[1,z[t-1]]&& k[1,z[t]]==k[1,z[t-2]])//AR(2)
{lp_p1 =lp + normal_lpdf(y[1,t] |(alphas[1,z[t]]+betas1[1,z[t]]*y[1,t-1]+betas2[1,z[t]]*y[1,t-2]), sigmas[1,z[t]]);}

else if (k[1,z[t]]==3 && k[1,z[t]]==k[1,z[t-1]]&& k[1,z[t]]==k[1,z[t-2]]&& k[1,z[t]]==k[1,z[t-3]])//AR(3)
{lp_p1 =lp + normal_lpdf(y[1,t] |(alphas[1,z[t]]+betas1[1,z[t]]*y[1,t-1]+betas2[1,z[t]]*y[1,t-2]+betas3[1,z[t]]*y[1,t-3]), sigmas[1,z[t]]);}
lp = lp_p1;

//Acy
if (k[2,z[t]]==1 && k[2,z[t]]==k[2,z[t-1]])//AR(1)
  {lp_p1 =lp + normal_lpdf(y[2,t] |(alphas[2,z[t]]+betas1[2,z[t]]*y[2,t-1]), sigmas[2,z[t]]);
  }
else if (k[2,z[t]]==2 && k[2,z[t]]==k[2,z[t-1]]&& k[2,z[t]]==k[2,z[t-2]])//AR(2)
{lp_p1 =lp + normal_lpdf(y[2,t] |(alphas[2,z[t]]+betas1[2,z[t]]*y[2,t-1]+betas2[2,z[t]]*y[2,t-2]), sigmas[1,z[t]]);}

else if (k[2,z[t]]==3 && k[2,z[t]]==k[2,z[t-1]]&& k[2,z[t]]==k[2,z[t-2]]&& k[2,z[t]]==k[2,z[t-3]])//AR(3)
{lp_p1 =lp + normal_lpdf(y[2,t] |(alphas[2,z[t]]+betas1[2,z[t]]*y[2,t-1]+betas2[2,z[t]]*y[2,t-2]+betas3[2,z[t]]*y[2,t-3]), sigmas[2,z[t]]);}
lp = lp_p1;


//Acz
if (k[3,z[t]]==1 && k[3,z[t]]==k[3,z[t-1]])//AR(1)
  {lp_p1 =lp + normal_lpdf(y[3,t] |(alphas[3,z[t]]+betas1[3,z[t]]*y[3,t-1]), sigmas[3,z[t]]);
  }
else if (k[3,z[t]]==2 && k[3,z[t]]==k[3,z[t-1]]&& k[3,z[t]]==k[3,z[t-2]])//AR(2)
{lp_p1 =lp + normal_lpdf(y[3,t] |(alphas[3,z[t]]+betas1[3,z[t]]*y[3,t-1]+betas2[3,z[t]]*y[3,t-2]), sigmas[3,z[t]]);}

else if (k[3,z[t]]==3 && k[3,z[t]]==k[3,z[t-1]]&& k[3,z[t]]==k[3,z[t-2]]&& k[3,z[t]]==k[3,z[t-3]])//AR(3)
{lp_p1 =lp + normal_lpdf(y[3,t] |(alphas[3,z[t]]+betas1[3,z[t]]*y[3,t-1]+betas2[3,z[t]]*y[3,t-2]+betas3[3,z[t]]*y[3,t-3]), sigmas[3,z[t]]);}
lp = lp_p1;
}

for (m in NAS[2:N]) { // looping over NAS
  lp_p1 =lp +  log_theta_tr[z[m-1],z[m]] + poisson_lpmf(u[count]|lambda[z[m]]);
  lp = lp_p1;
  count=count+1;
}

target += lp;
}
