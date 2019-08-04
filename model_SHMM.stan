data {
  int<lower=1> K;  // num categories
  int<lower=0> T;  // num instances
  int<lower=0> N;   // num of NAS
  int<lower=1,upper=K> z[T]; // behaviours
  int<lower=0> u[T]; // sojuorn times
  int NAS[N]; // indexes of NAS points
  matrix[3,T] y;// obs data
  //real y[T];// obs data
  int <lower=0,upper=3> k[3*K, 4];//AR() orders
 //int <lower=0,upper=2> kz[K]; //ACCX AR() orders
  }
parameters {
  //sojourn times parameters
  real <lower=0> lambda[K];  
  //autorregresive models.
  real betas1[3,K];
  real betas2[3,K];
  real betas3[3,K];
  real alphas[3,K];
  real<lower=0> sigmas[3,K];
  // K x K tpm
  simplex[K] theta[K]; 
}


transformed parameters{
matrix[K, K] ta; 
for(j in 1:K){
for(i in 1:K){
ta[i,j]= theta[i,j];
}
}
for(i in 1:K){
ta[i,i]= 0;
}
}

model {
  vector[K] log_theta_tr[K];
  real lp;
  real lp_p1;
  real lp_px;
  real lp_py;
  real lp_pz;
      
  
  // prior for alphas and betas ??
  //betas1 ~ student_t(3, 0, 1);
  //betas2 ~ student_t(3, 0, 1);
  //alphas ~ student_t(3, 0, 1);
  for(i in 1:K)
  {lambda[i]~normal(15, 5);
    for(j in 1:3)
      sigmas[j,i] ~ cauchy(0,2);
  }
// transpose the tpm and take natural log of entries
for (n_from in 1:K)
  for (n in 1:K)
  if(n_from !=n)
    log_theta_tr[n, n_from] = log(ta[ n,n_from]);
// Compute CDL

//First obs
lp = log(0.5) + poisson_lpmf(u[1]|lambda[z[1]])+normal_lpdf(y[1,1] | 0, sigmas[1,z[1]])+
normal_lpdf(y[2,1] | 0, sigmas[2,z[1]])+normal_lpdf(y[3,1] | 0, sigmas[3,z[1]]);

for (t in 4:T) { // looping over all observations pdfs
//Acc x y z
lp_px = normal_lpdf(y[1,t] |(alphas[1,z[t]]+betas1[1,z[t]]*y[1,t-1]+betas2[1,z[t]]*y[1,t-2]*k[z[t],3]+betas3[1,z[t]]*y[1,t-3]*k[z[t],4]), sigmas[1,z[t]]);

lp_py = normal_lpdf(y[2,t] |(alphas[2,z[t]]+betas1[2,z[t]]*y[2,t-1]+betas2[2,z[t]]*y[2,t-2]*k[5+z[t],3]+betas3[2,z[t]]*y[2,t-3]*k[5+z[t],4]), sigmas[2,z[t]]);

lp_pz = normal_lpdf(y[3,t] |(alphas[3,z[t]]+betas1[3,z[t]]*y[3,t-1]+betas2[3,z[t]]*y[3,t-2]*k[10+z[t],3]+betas3[3,z[t]]*y[3,t-3]*k[10+z[t],4]), sigmas[3,z[t]]);

lp = lp+lp_px+lp_py+lp_pz;

}

for (m in NAS[2:N]) { // looping over NAS
 lp_p1 = log_theta_tr[z[m-1],z[m]]+ poisson_lpmf(u[m]|lambda[z[m]]);
 lp = lp+lp_p1;
  //count=count+1;
}

target += lp;
}
