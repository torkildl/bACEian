data {

	int <lower=0> n_obs;    // number of observations
	int <lower=0> n_fam;    // number of families
	int <lower=0> n_famtw;  // number of nested twins in families
	
	real y[n_obs];    // responses
	int mz[n_obs];    // mz dummy
	int dz[n_obs];    // dz dummy
	int fam[n_obs];   // family indicator (1:nfam)
	int famtw[n_obs]; // nested twin indicator (1:nfamtw)
	
}

transformed data {
  
  real z1[n_obs];
  real z2[n_obs];
  
  // compute weights
  for(i in 1:n_obs) {
      z1[i] = sqrt(.5)*dz[i];
      z2[i] = mz[i] + sqrt(.5)*dz[i];
  }
  
}

parameters {
	
	// mean
	real mu;
	
	// "random-effects" sd
	real<lower=0> a; 
	real<lower=0> c; 
	real<lower=0> e; 
	
	// "random effects"
	vector[n_famtw] a_unique;
	vector[n_fam] a_shared;
	vector[n_fam] c_shared;
	
}

transformed parameters {

	vector[n_obs] yhat;

  // predicted responses
	for(i in 1:n_obs) {
		yhat[i] = mu + a_unique[ famtw[i] ]*z1[i] + a_shared[ fam[i] ]*z2[i] + c_shared[ fam[i] ];
	}
	
}

model {

	// prior distributions
	mu ~ normal(0, 5);
	a_unique ~ normal(0, a);
	a_shared ~ normal(0, a);
	c_shared ~ normal(0, c);

	// likelihood
	y ~ normal(yhat, e);
	
}

generated quantities {

	real A;
	real C;
	real E;
	real V;
	real Asd;
	real Csd;
	real Esd;

	// variances
	A = a^2; 
	C = c^2;
	E = e^2;
	V = A + C + E;
	
	// standardized variance components
	Asd = A / V; // narrow-sense heritability
	Csd = C / V;
	Esd = E / V;

}
