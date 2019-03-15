// Stan code for mixed effects-type ACE model with Dirichlet prior
//
// 

data {
  int <lower=0> n_fam; // number of families
  int <lower=0> n_famtw_mz; // number of mz nested twins in families
  int <lower=0> n_famtw_dz; // number of dz nested twins in families
  real y_mz[n_famtw_mz * 2]; // responses
  real y_dz[n_famtw_dz * 2]; // responses
  int fam_mz[n_famtw_mz * 2]; // family indicator (1:nfam)
  int fam_dz[n_famtw_dz * 2]; // family indicator (1:nfam)
  real<lower = 0> outcome_sd;
  real outcome_mean;
}

transformed data {
  vector[3] dirichlet_prior;
  dirichlet_prior = rep_vector(1, 3);
}

parameters {
  // mean
  real mu;
  // overall variance
  real<lower = 0> sigma;
  // variance components
  simplex[3] var_comp_shares;
  // "random-effects" sd for genetics
  vector[n_fam] a_shared_std;
  // random effects sd for common env
  vector[n_fam] c_shared_std;
}

transformed parameters {
  // "random effects"
  vector[n_fam] a_shared;
  vector[n_fam] c_shared;
  real A;
  real C;
  real E;
  // real V;
  real Asd;
  real Csd;
  real Esd;
  real a;
  real c;
  real e_sigma;
  Asd = var_comp_shares[1];
  Csd = var_comp_shares[2];
  Esd = var_comp_shares[3];
  A = Asd * sigma^2;
  C = Csd * sigma^2;
  E = Esd * sigma^2;
  a = sqrt(A);
  c = sqrt(C);
  e_sigma = sqrt(E);
  a_shared = a * a_shared_std;
  c_shared = c * c_shared_std;
}

model {
  vector[n_famtw_mz * 2] y_mz_expected;
  vector[n_famtw_dz * 2] y_dz_expected;
  mu ~ normal(outcome_mean, outcome_mean * 0.2);
  sigma ~ normal(outcome_sd, outcome_sd * 0.3);
  // a ~ normal(0, 5);
  // c ~ normal(0, 5);
  // e_sigma ~ normal(0, 5);
  var_comp_shares ~ dirichlet(dirichlet_prior);
  a_shared_std ~ normal(0,1);
  c_shared_std ~ normal(0,1);
  // model
  for (i in 1:(n_famtw_mz * 2)){
    y_mz_expected[i] = mu + a_shared[fam_mz[i]] + c_shared[fam_mz[i]];
  }
  for (i in 1:(n_famtw_dz * 2)){
    y_dz_expected[i] = mu + sqrt(0.5) * a_shared[fam_dz[i]] + c_shared[fam_dz[i]];
  }
  target += normal_lpdf(y_mz | y_mz_expected, e_sigma);
  target += normal_lpdf(y_dz | y_dz_expected, sqrt(E + 0.5 * A));
}

generated quantities {
  vector[n_famtw_mz * 2] y_rep_mz;
  vector[n_famtw_dz * 2] y_rep_dz;
  for (i in 1:(n_famtw_mz * 2)){
    y_rep_mz[i] = normal_rng(mu + a_shared[fam_mz[i]] + c_shared[fam_mz[i]], e_sigma);
  }
  for (i in 1:(n_famtw_dz * 2)){
    y_rep_dz[i] = normal_rng( mu + sqrt(0.5) * a_shared[fam_dz[i]] + c_shared[fam_dz[i]], sqrt(E + 0.5 * A));
  }
}
