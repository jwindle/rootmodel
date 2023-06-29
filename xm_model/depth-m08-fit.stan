/* -*- c-basic-offset: 4; -*- */
data {
	int<lower=1> N;
	int<lower=1> M;
	real<lower=0> r;
	real <upper=0.0> min_depth;
	real <upper=0.0> max_depth;
	real mu_y_prior_mean;
	real mu_y_prior_std;
	real sig_y_prior_mean;
	real sig_y_prior_std;
	int<lower=1, upper=M> group[N]; /* group membership */
	real y[N];
	/* /\* Below is only necessary for going back to mu_dm and sig_dm *\/ */
	/* real <lower=0> mu_dx; */
	/* int<lower=1> K; /\* ceil(r / mu_dx) - 1; does not work in stan *\/ */
}
transformed data {
	real lb = min_depth;
	real ub = max_depth;
	/* /\* Below is only necessary for going back to mu_dm and sig_dm *\/ */
	/* vector[K+1] rvec = rep_vector(r, K+1); */
	/* vector[K+1] dx = rep_vector(mu_dx, K+1); */
	/* dx[1] = 0.0; */
	/* vector[K+1] x = cumulative_sum(dx); */
	/* vector[K+1] x_left = fdim(rvec, x); */
	/* real x_left_one = sum(x_left); */
	/* real x_left_norm = sqrt(sum(square(x_left))); */
}
parameters {
	real mu_y[M];
	real<lower=0> sig_y[M];
}
model {
	for (i in 1:M) {
		mu_y[i] ~ normal(mu_y_prior_mean, mu_y_prior_std);
		sig_y[i] ~ normal(sig_y_prior_mean, sig_y_prior_std) T[0,];
	}
	for (i in 1:N) {
		int g_i = group[i];
		y[i] ~ normal(mu_y[g_i], sig_y[g_i]) T[lb, ub];
	}
}
/* generated quantities { */
/* 	real mu_dm[M]; */
/* 	real sig_dm[M]; */
/* 	for (i in 1:M) { */
/* 		mu_dm[i] = mu_y[i] / x_left_one; */
/* 		sig_dm[i] = sig_y[i] / x_left_norm; */
/* 	} */
/* } */
