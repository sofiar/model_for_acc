data {
  int<lower=0> N;
  int<lower=0> M;
  vector[N] y;
  int ydims[M+1];
  int rdims[M+1];
  
  
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}

model {
  sigma ~ student_t(3, 0, 1);
  beta~normal(0, 10);//out?
  for(i in 1:M)
  y[(ydims[i] + 1):(ydims[i+1]-1)] ~ normal(alpha + beta * y[ydims[i]:(ydims[i+1] - 2)], sigma);
  }

generated quantities {
vector[(N-M)] rss;
vector[N] ypred;
  ypred[1] = y[1];
  for(i in 1:M)
  {
  rss[(rdims[i]+1):(rdims[i+1]) ]= y[(ydims[i] + 1):(ydims[i+1]-1)]-(alpha + beta * y[ydims[i]:(ydims[i+1] - 2)]);
  }
  
  for (n in 2:N)
  {
  ypred[n]=normal_rng(alpha + beta * ypred[n-1], sigma);
  }
  }
  