data {
  int<lower=1> K;  // num categories
  int<lower=0> T;  // num instances
  int<lower=0> N;   // num of NAS
  int<lower=1,upper=K> z[T]; // behaviours
  int<lower=0> u[T]; // sojuorn times
  int NAS[N]; // indexes of NAS points
  matrix[3,T] y;// obs data
  
  //int <lower=0,upper=3> k[3*K, 4];//AR() orders
 //int <lower=0,upper=2> kz[K]; //ACCX AR() orders
  }
parameters {
  //sojourn times parameters
  real <lower=0> lambda[K];  
  //autorregresive models.
  real alphas[3,K];
  real betas1[3,K];
  real betas2[6]; //fixed from the sysmtem
  real betas3[2]; //fixed from the sysmtem 
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
      
  
  // priors 
  for (i in 1:3)
  {
  alphas[i] ~ student_t(10, 0, 1);  
  betas1[i] ~ student_t(3, 0, 1);
  }
  
  //betas2 ~ student_t(3, 0, 1);
  
  
  
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

if(z[t]==1) // Baha 1
{
//Acc x --> AR(1)
lp_px = normal_lpdf(y[1,t] |alphas[1,z[t]]+betas1[1,z[t]]*y[1,t-1], sigmas[1,z[t]]);
//Acc y --> AR(2)
lp_py =  normal_lpdf(y[2,t] |alphas[2,z[t]]+betas1[2,z[t]]*y[2,t-1]+betas2[1]*y[2,t-2], sigmas[2,z[t]]); 
//Acc z --> AR(1)
lp_pz =  normal_lpdf(y[3,t] |alphas[3,z[t]]+betas1[3,z[t]]*y[3,t-1], sigmas[3,z[t]]); 
}

if(z[t]==2)// Beha 2
{
//Acc x --> AR(1)
lp_px = normal_lpdf(y[1,t] |alphas[1,z[t]]+betas1[1,z[t]]*y[1,t-1], sigmas[1,z[t]]);
//Acc y --> AR(2)
lp_py =  normal_lpdf(y[2,t] |alphas[2,z[t]]+betas1[2,z[t]]*y[2,t-1]+betas2[2]*y[2,t-2], sigmas[2,z[t]]); 
//Acc z --> AR(3)
lp_pz =  normal_lpdf(y[3,t] |alphas[3,z[t]]+betas1[3,z[t]]*y[3,t-1]+betas2[3]*y[3,t-2]+betas3[1]*y[3,t-3], sigmas[3,z[t]]); 
}

if(z[t]==3)// Beha 3
{
//Acc x --> AR(1)
lp_px = normal_lpdf(y[1,t] |alphas[1,z[t]]+betas1[1,z[t]]*y[1,t-1], sigmas[1,z[t]]);
//Acc y --> AR(1)
lp_py =  normal_lpdf(y[2,t] |alphas[2,z[t]]+betas1[2,z[t]]*y[2,t-1],sigmas[2,z[t]]); 
//Acc z --> AR(1)
lp_pz =  normal_lpdf(y[3,t] |alphas[3,z[t]]+betas1[3,z[t]]*y[3,t-1],sigmas[3,z[t]]); 
}

if(z[t]==4)// Beha 4
{
//Acc x --> AR(2)
lp_px = normal_lpdf(y[1,t]|alphas[1,z[t]]+betas1[1,z[t]]*y[1,t-1]+betas2[4]*y[1,t-2], sigmas[1,z[t]]);
//Acc y --> AR(3)
lp_py =  normal_lpdf(y[2,t]|alphas[2,z[t]]+betas1[2,z[t]]*y[2,t-1]+betas2[5]*y[2,t-2]+betas3[2]*y[2,t-3],sigmas[2,z[t]]); 
//Acc z --> AR(1)
lp_pz =  normal_lpdf(y[3,t]|alphas[3,z[t]]+betas1[3,z[t]]*y[3,t-1],sigmas[3,z[t]]); 
}

if(z[t]==5)// Beha 5
{
//Acc x --> AR(2)
lp_px = normal_lpdf(y[1,t]|alphas[1,z[t]]+betas1[1,z[t]]*y[1,t-1]+betas2[6]*y[1,t-2], sigmas[1,z[t]]);
//Acc y --> AR(1)
lp_py =  normal_lpdf(y[2,t]|alphas[2,z[t]]+betas1[2,z[t]]*y[2,t-1],sigmas[2,z[t]]); 
//Acc z --> AR(1)
lp_pz =  normal_lpdf(y[3,t]|alphas[3,z[t]]+betas1[3,z[t]]*y[3,t-1],sigmas[3,z[t]]); 
}


lp = lp+lp_px+lp_py+lp_pz;

}

for (m in NAS[2:N]) { // looping over NAS
 lp_p1 = log_theta_tr[z[m-1],z[m]]+ poisson_lpmf(u[m]|lambda[z[m]]);
 lp = lp+lp_p1;
  //count=count+1;
}

target += lp;
}
