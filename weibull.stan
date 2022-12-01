functions{
  
  vector path(int n, vector hours, int[] units, vector theta){
    vector[n] D;
    for(i in 1:n){
      D[i] = hours[i]/theta[units[i]];
    }
    return(D);
  }
}

data {
  int<lower=1> L;
  int<lower=1> n;
  real<lower=0> C;
  int<lower=1> units[n];
  vector[n] y;
  vector[n] hours;
  real<lower=0> Df;
  real<lower=0, upper = 1> alpha;

  real<lower=0> a_lambda;
  real<lower=0> b_lambda;
  real<lower=0> a_tau;
  real<lower=0> b_tau;
}


parameters {
  vector<lower=0>[L] theta;
  real<lower=0> beta;
  real<lower=0> lambda;
  real<lower=0> tau;
}

transformed parameters{
  real<lower = 0> sigma = sqrt(1/tau);
}


model{
  vector[n] D = path(n, hours, units, theta);
  y ~ normal(D, sigma);
  theta ~ weibull(beta, pow(lambda, -1/beta) );
  lambda ~ gamma(a_lambda, b_lambda);
  tau ~ gamma(a_tau, b_tau);
}

generated quantities{
  real<lower = 0> r4500 = exp(-(lambda/pow(Df, beta))*pow(4500.0/C, beta));
  real<lower = 0> talpha = pow(( -(pow(Df, beta)/lambda)*log(1-alpha)), (1/beta));
}

