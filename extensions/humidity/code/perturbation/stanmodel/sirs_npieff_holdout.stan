// Paired Baseline (K = 0) / humidity_seasonal (K = 3).
// Same SIRS equations and common priors as stanmodel/sirs_npi.stan;
// npi_anchor = 0 reproduces the Shikoku prior exception in sirs_npi_tweak.stan.
// Only training cases and training covariates enter the model.
// N_total fixes the latent-path prior horizon: here 2013--2019, NOT 2013--2023.
// Future npi values are sampled jointly under the original finite-path prior,
// with no holdout cases and no holdout SIRS states evaluated in this model.

data {
  int<lower=2> N_train;
  int<lower=N_train> N_total;
  int<lower=0, upper=N_total> Nnonpi;
  int<lower=0, upper=1> npi_anchor;
  int<lower=0> K;
  array[N_train] int<lower=1, upper=52> week;
  array[N_train] int<lower=0> cases;
  matrix[N_train, K] X_train;
  vector<lower=0>[K] theta_sd;
  real<lower=0> mu;
  real<lower=0> pop;
  real<lower=0> gamma;
}

transformed data {
  int<lower=0> M = N_total - Nnonpi;
  int<lower=0, upper=1> has_npi = M > 0;
  int<lower=0> N_increments = M > 0 ? M - 1 : 0;
}

parameters {
  real<lower=0, upper=1> S0;
  real<lower=0, upper=1> I0;
  real<lower=0, upper=1> rho;
  array[52] real<lower=0> beta;
  real<lower=0> sigma;
  real<lower=0> phi;
  real<lower=0> tau;
  real<lower=0> omega;
  real<lower=0> sigma_npi;
  array[has_npi] real<lower=0> npi_first;
  vector[N_increments] z_npi;
  vector[K] theta;
}

transformed parameters {
  vector<lower=0>[N_train] S;
  vector<lower=0>[N_train] I;
  vector<lower=0>[N_train] R;
  vector<lower=0>[N_train] C;
  real<lower=0> delta = 1 / tau;
  vector<lower=0>[M] npi;
  vector<lower=0>[N_total] npieff = rep_vector(1, N_total);

  // Exact noncentered version of the original positive npi random walk.
  // The transformed-parameter constraint retains the original positive support.
  if (M > 0) {
    npi[1] = npi_first[1];
    if (M > 1) {
      for (j in 2:M) {
        npi[j] = npi[j - 1] + sigma_npi * z_npi[j - 1];
      }
    }
    for (j in 1:M) {
      npieff[Nnonpi + j] = npi[j];
    }
  }

  S[1] = S0 * pop;
  I[1] = I0 * pop;
  R[1] = (1 - S0 - I0) * pop;
  C[1] = 0;

  for (i in 2:N_train) {
    real humidity_multiplier = 1;
    real foi;
    real Sout;
    real StoI;
    real Iout;
    real ItoR;
    real Rout;
    real RtoS;
    if (K > 0) {
      humidity_multiplier = exp(X_train[i] * theta);
    }
    foi = beta[week[i]] * humidity_multiplier
          * (I[i - 1] + omega) / pop * npieff[i];
    Sout = (1 - exp(-(foi + mu))) * S[i - 1];
    StoI = foi / (foi + mu) * Sout;
    Iout = (1 - exp(-(gamma + mu))) * I[i - 1];
    ItoR = gamma / (gamma + mu) * Iout;
    Rout = (1 - exp(-(delta + mu))) * R[i - 1];
    RtoS = delta / (delta + mu) * Rout;

    S[i] = S[i - 1] + RtoS - Sout + mu * pop;
    I[i] = I[i - 1] + StoI - Iout;
    R[i] = R[i - 1] + ItoR - Rout;
    C[i] = StoI * rho;
  }
}

model {
  beta ~ normal(8, 2);
  for (i in 2:52) {
    beta[i] ~ normal(beta[i - 1], sigma);
  }
  beta[1] ~ normal(beta[52], sigma);

  // Under npi[j] = npi[j-1] + sigma_npi*z_npi[j-1],
  // |d npi[2:M] / d z_npi| = sigma_npi^(M-1).
  // This cancels the transition normals' sigma_npi^(-(M-1)) factors.
  // Do not add a second random-walk density or an additional Jacobian here.
  if (M > 1) {
    z_npi ~ std_normal();
  }
  if (M > 0 && npi_anchor == 1) {
    npi ~ normal(1, 0.2);
  }
  sigma ~ normal(0, 1);
  sigma_npi ~ normal(0, 0.1);
  tau ~ normal(0, 200);
  S0 ~ normal(0.08, 0.02);
  I0 ~ normal(0, 0.001);
  rho ~ normal(0, 0.02);
  phi ~ normal(0, 10);
  omega ~ normal(0, 200);

  if (K > 0) {
    theta ~ normal(0, theta_sd);
  }
  for (i in 2:N_train) {
    cases[i] ~ neg_binomial_2(C[i], phi);
  }
}
