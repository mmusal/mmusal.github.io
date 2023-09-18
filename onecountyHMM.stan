data { 
  int<lower=1> T;//number of time periods
  int K;
  int START;
  int END;
  int<lower=0> N;//number of spatial areas
  int<lower=0>  y[T,N];              // count outcomes
  matrix [T,N] log_E;
}
transformed data {
  int ITER=END-START+1;
}
parameters {
ordered[K] mu[1]; 
  //real<lower=0> sigma[K];
  simplex[K] pi1[1];
  simplex[K] A[K,1];
//  real <lower=0> mu_sd;
//  real mu_mean;
  real epsilon[T];
}
transformed parameters {
matrix[K,1] lambda [T] ;
real logalpha[T,K];
real accumulator[K];
    for(t in 1:T){
  for(k in 1:K){     
      lambda[t][k,1]=exp(
        log_E[t,1]+mu[1,k]+epsilon[t]); 
                }
                        }
for(k in 1:K){

logalpha[1,k] = log(pi1[1,k]) + poisson_lpmf(y[1,1] | lambda[1,k,1]);

}

for (t in 2:T) {
for (j in 1:K) { // j = current (t)

accumulator[1] = logalpha[t-1, 1] + log(A[1, 1,j]) + poisson_lpmf(y[t,1] | lambda[t,j,1]);
accumulator[2] = logalpha[t-1, 2] + log(A[2, 1,j]) + poisson_lpmf(y[t,1] | lambda[t,j,1]);

//I do not like this
logalpha[t, j] = log_sum_exp(to_vector(accumulator));
}
}
} // Forward

model {
 target += log_sum_exp(logalpha[T]);
 //this is where the issue is

  mu[1]~normal(0,3);


epsilon~normal(0,3);
}

