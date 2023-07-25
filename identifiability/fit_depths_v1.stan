/* -*- c-basic-offset: 4; -*- */
/* Jesse Windle, 2023 */
data {
	int<lower=1> N;              /* number of total observations       */
	int<lower=1> M;              /* number of groups                   */
	int<lower=1> K;              /* number of kinks                    */
	int<lower=1, upper=M> group[N]; /* group membership */
	real<lower=0> r;             /* radius at which a root is detected */
	real mu_dx_prior_mean;
	real mu_dx_prior_std;
	real mu_dm_prior_mean;
	real mu_dm_prior_cond_slope;
	real mu_dm_prior_cond_std;
	real sig_dx_prior_mean;
	real sig_dx_prior_std;
	real sig_dm_prior_mean;
	real sig_dm_prior_std;
	real y[N];
}
transformed data {
	vector[K] rvec = rep_vector(r, K);
}
parameters {
	real mu_dm[M];
	real<lower=0> sig_dm[M];
	real<lower=0> mu_dx[M];
	real<lower=0> sig_dx[M];
	real dm_0[N];
	vector[K] dm[N];
	/* vector<lower=0>[K] dx[N]; */
	vector<lower=r/(K+1.5), upper=r/(1.5)>[K] dx[N];
}
transformed parameters {
	real mu_dm_cond_mean[M];
	vector[K] x[N];
	vector[K] x_left[N];
	real z[N];
	for (i in 1:M) {
		mu_dm_cond_mean[i] = mu_dm_prior_mean + mu_dm_prior_cond_slope * (mu_dx[i] - mu_dx_prior_mean);
	}
	for (i in 1:N) {
		x[i] = cumulative_sum(dx[i]);
		x_left[i] = fdim(rvec, x[i]);
		z[i] = dm_0[i] * r + sum(dm[i] .* x_left[i]);
	}
}
model {
	for (i in 1:M) {
		mu_dx[i] ~ normal(mu_dx_prior_mean, mu_dx_prior_std) T[0,];
		mu_dm[i] ~ normal(mu_dm_cond_mean, mu_dm_prior_cond_std);
		sig_dx[i] ~ normal(sig_dx_prior_mean, sig_dx_prior_std) T[0, r];
		sig_dm[i] ~ normal(sig_dm_prior_mean, sig_dm_prior_std) T[0,];
	}
	for (i in 1:N) {
		int j = group[i];
		dm_0[i] ~ normal(mu_dm[j], sig_dm[j]);
		dm[i] ~ normal(mu_dm[j], sig_dm[j]);
		for (k in 1:K) {
			/* dx[i][k] ~ normal(mu_dx[j], sig_dx[j]) T[0,]; */
			dx[i][k] ~ normal(mu_dx[j], sig_dx[j]) T[r/(K+1.5),r/(1.5)];
		}
		y[i] ~ normal(z[i], 0.01);
	}
}
