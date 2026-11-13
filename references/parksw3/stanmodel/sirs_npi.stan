data {
  int N;
  int Nnonpi;
  int week[N];
  int cases[N];
  real mu;
  real pop;
  real gamma;
}

parameters {
  real<lower=0, upper=1> S0;
  real<lower=0, upper=1> I0;
  real<lower=0, upper=1> rho;
  real<lower=0> beta[52];
  real<lower=0> sigma;
  real<lower=0> phi;
  real<lower=0> tau;
  real<lower=0> omega;
  real<lower=0> npi[N-Nnonpi];
  real<lower=0> sigma_npi;
}

transformed parameters {
  vector<lower=0>[N] S;
  vector<lower=0>[N] I;
  vector<lower=0>[N] R;
  vector<lower=0>[N] C;
  real<lower=0> delta = 1/tau;
  vector<lower=0>[N] npieff;

  S[1] = S0 * pop;
  I[1] = I0 * pop;
  R[1] = (1-S0-I0) * pop;
  C[1] = 0;

  for (i in 1:Nnonpi) {
    npieff[i] = 1;
  }

  for (i in (Nnonpi+1):N) {
    npieff[i] = npi[i-Nnonpi];
  }

  for (i in 2:N) {
    real foi = beta[week[i]] * (I[i-1]+omega)/pop * npieff[i];
    real Sout = (1-exp(-(foi+mu))) * S[i-1];
    real StoI = foi/(foi+mu) * Sout;
    real Iout = (1-exp(-(gamma+mu))) * I[i-1];
    real ItoR = gamma/(gamma+mu) * Iout;
    real Rout = (1-exp(-(delta+mu))) * R[i-1];
    real RtoS = delta/(delta+mu) * Rout;

    S[i] = S[i-1] + RtoS - Sout + mu * pop;
    I[i] = I[i-1] + StoI - Iout;
    R[i] = R[i-1] + ItoR - Rout;
    C[i] = StoI * rho;
  }
}

model {
  beta ~ normal(8, 2);
  for (i in 2:52) {
    beta[i] ~ normal(beta[i-1], sigma);
  }
  beta[1] ~ normal(beta[52], sigma);

  npi ~ normal(1, 0.2);
  for (i in 2:(N-Nnonpi)) {
    npi[i] ~ normal(npi[i-1], sigma_npi);
  }

  sigma ~ normal(0, 1);
  sigma_npi ~ normal(0, 0.1);

  tau ~ normal(0, 200);

  S0 ~ normal(0.08, 0.02);
  I0 ~ normal(0, 0.001);
  rho ~ normal(0, 0.02);
  phi ~ normal(0, 10);
  omega ~ normal(0, 200);

  for (i in 2:N) {
    cases[i] ~ neg_binomial_2(C[i], phi);
  }
}
