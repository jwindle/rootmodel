/* -*- c-basic-offset: 4; -*- */
data {
	int<lower=1> N;              /* number of total observations       */
	int<lower=1> M;              /* number of groups                   */
	real<lower=0> mu_dx;         /* dx step size */
	int<lower=1> K;              /* number of kinks                    */
	int<lower=1, upper=M> group[N]; /* group membership */
	real<lower=0> r;             /* radius at which a root is detected */
	real mu_dm_prior_mean;
	real mu_dm_prior_std;
	real sig_dm_prior_mean;
	real sig_dm_prior_std;
	real shape_dm_prior_mean;
	real shape_dm_prior_std;
	real y[N];
}
transformed data {
	vector[K] rvec = rep_vector(r, K);
	vector[K] dx = rep_vector(mu_dx, K);
	vector[K] x = cumulative_sum(dx);
	vector[K] x_left = fdim(rvec, x);
	real sig_dm_alpha = square(sig_dm_prior_mean / sig_dm_prior_std) + 2;
	real sig_dm_beta = sig_dm_prior_mean * (sig_dm_alpha - 1);

}
parameters {
	real mu_dm[M];
	real<lower=0> sig_dm[M];
	real shape_dm[M];
	real dm_0[N];
	vector[K] dm[N];
}
transformed parameters {
	real z[N];
	for (i in 1:N) {
		z[i] = dm_0[i] * r + sum(dm[i] .* x_left);
	}
}
model {
	for (i in 1:M) {
		mu_dm[i] ~ normal(mu_dm_prior_mean, mu_dm_prior_std);
		sig_dm[i] ~ normal(sig_dm_prior_mean, sig_dm_prior_std) T[0,];
		/* sig_dm[i] ~ inv_gamma(sig_dm_alpha, sig_dm_beta); */
		shape_dm[i] ~ normal(shape_dm_prior_mean, shape_dm_prior_std);
	}
	for (i in 1:N) {
		int j = group[i];
		dm_0[i] ~ skew_normal(mu_dm[j], sig_dm[j], shape_dm[j]);
		dm[i] ~ skew_normal(mu_dm[j], sig_dm[j], shape_dm[j]);
		y[i] ~ normal(z[i], 0.01);
	}
}
