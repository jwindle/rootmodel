/* -*- c-basic-offset: 4; -*- */
data {
	int<lower=1> N;
	int<lower=1> M;
	real<lower=0.> r;
	real<upper=0.> min_depth;
	real<upper=0.> max_depth;
	real <lower=0.> mu_dx;
	int<lower=1> K; /* ceil(r / mu_dx) - 1; does not work in stan */
	real mu_dm_prior_mean;
	real mu_dm_prior_std;
	real sig_dm_prior_mean;
	real sig_dm_prior_std;
	int<lower=1, upper=M> group[N]; /* group membership */
	real y[N];
}
transformed data {
	vector[K+1] rvec = rep_vector(r, K+1);
	vector[K+1] dx = rep_vector(mu_dx, K+1);
	dx[1] = 0.0;
	vector[K+1] x = cumulative_sum(dx);
	vector[K+1] x_left = fdim(rvec, x);
	real<upper=0.0> lb = min_depth / sum(x_left);
	real<upper=0.0> ub = max_depth / sum(x_left);
}
parameters {
	real mu_dm[M];
	real<lower=0> sig_dm[M];
	vector<lower=lb, upper=ub>[K+1] dm[N];
}
transformed parameters {
	real z[N];
	for(i in 1:N) {
		z[i] = dot_product(dm[i], x_left);
	}
}
model {
	for (i in 1:M) {
		mu_dm[i] ~ normal(mu_dm_prior_mean, mu_dm_prior_std);
		sig_dm[i] ~ normal(sig_dm_prior_mean, sig_dm_prior_std) T[0,];
		/* sig_dm[i] ~ inv_gamma(sig_dm_alpha, sig_dm_beta); */
	}
	for (i in 1:N) {
		int g_i = group[i];
		for (j in 1:(K+1)) {
			dm[i][j] ~ normal(mu_dm[g_i], sig_dm[g_i]) T[lb, ub];
		}
		y[i] ~ normal(z[i], 0.01);
	}
}
generated quantities {
	vector[K+1] m[N];
	vector[K+1] theta[N];
	vector[K+1] dl[N];
	for (i in 1:N) {
		m[i][1] = dm[i][1];
		for (j in 1:K) {
			m[i][j+1] = m[i][j] + dm[i][j+1];
		}
		theta[i] = atan(m[i]);
		dl[i] = dx ./ cos(theta[i]);
	}
}
