data {
  int<lower=1> K;  // num categories
  int<lower=0> T;  // num instances
  int<lower=0> N;  // num of NAS
  int<lower=1,upper=K> z[T]; // behaviours
  int<lower=0> u[N]; // sojuorn times
  int NAS[N]; // indexes of NAS points
  //matrix[T, 3] y;// obs data
  real y[T];// obs data
  int <lower=0,upper=2> k[K]; //AR() orders
  }
parameters {
  real <lower=0> lambda[K]; //sojourn times parameters 
  //autorregresive models. Assume AR(1) for now
  real betas1[K];
  real betas2[K];// automatic (?)
  real alphas[K];
  real<lower=0> sigmas[K];
  simplex[K] theta[K]; // K x K tpm
}

model {
  vector[K] log_theta_tr[K];
  real lp;
  real lp_p1;
  int count=1;
  int k1=1;
  int k2=2;
  
  // prior for alphas and betas 
  betas1 ~ student_t(3, 0, 1);
  betas2 ~ student_t(3, 0, 1);
  alphas ~ student_t(3, 0, 1);
  sigmas ~ cauchy(0,2);
  
  
// transpose the tpm and take natural log of entries
for (n_from in 1:K)
  for (n in 1:K)
    log_theta_tr[n, n_from] = log(theta[n_from, n]);
// Compute CDL

//First obs
lp = log(0.5) + poisson_lpmf(u[1]|lambda[z[1]])+normal_lpdf(y[1] | 0, sigmas[z[1]]);


// ACA HAY QUER AGREGAR UN IF QUE TENGA QUE VER CON SI ES LA PRIMERA OBS DEL ESTADO
for (t in 3:T) { // looping over all observations pdfs
if (k[z[t]]==1)//AR(1)
{lp_p1 =lp + normal_lpdf(y[t] |(alphas[z[t]]+betas1[z[t]]*y[t-1]), sigmas[z[t]]);}
  
else if (k[z[t]]==2)//AR(2)
{lp_p1 =lp + normal_lpdf(y[t] |(alphas[z[t]]+betas1[z[t]]*y[t-1]+betas2[z[t]]*y[t-2]), sigmas[z[t]]);}

lp = lp_p1;
}

for (m in NAS[2:N]) { // looping over NAS
  lp_p1 =lp +  log_theta_tr[z[m-1],z[m]] + poisson_lpmf(u[count]|lambda[z[m]]);
  lp = lp_p1;
count=count+1;
}

target += lp;
}
