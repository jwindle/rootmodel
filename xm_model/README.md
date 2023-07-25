# Overview

Here we build models of root growth.  We model the trajectory of roots
given depth data, patterns in root growth emergence, and both of these
together.


# Root trajectories

Every model has two parts, one to simulate pseudo data and another to
fit the model.  Every pair corresponds to files named like
`depth-mXX-sim.stan` and `depth-mXX-fit.stan` where `XX` could be
replaced by two numbers, like `01`, `02`, etc.

These files are used by other scripts to compile and run models.
There are three different routines we can: 

1. simulate from a model (`depth-MAIN-sim_and_plot.R`),
2. fit simulated data (`depth-MAIN-fit_sim.R`), and
3. and fit real data.  

Fitting real data is broken into four parts:

1. Load data (`depth-MAIN-fit_real-laod.R`)
1. Fit data (`depth-MAIN-fit_real-fit.R`)
1. Generate posterior predictive distributions (`depth-MAIN-fit_real-pp.R`)
1. Compare models (`depth-MAIN-fit_real-compare.R`)

Each of these options is accompanied by a config file for controlling
which models are run and model set up.

The models are:

- 01: multivariate normal truncated to got go below a certain depth, no latents
- 02: multivariate normal truncated to got go below a certain depth, with latents
- 03: skew normal under standard parameterization
- 04: skew normal with centered parameterization, via skew norm latents
- 05: skew normal with centered parameterization, via parameters
- 06: truncated normal
- 08: like 01, but using a prior that is meaningful in terms of observed depth
- 09: normal with correlated changing slopes
- 10: skew normal
- 21: modeling angle, not changing slopes, using a logistic transform
- 22: modeling angle, using atan transform (and changing secant approximation)
- 23: like m10, but with dynamic parameters that come from a Gaussian process


# Root emergence

We model root emergence using via `fit_times.stan` which is run in
`fit_times.R`


# Simulating plant growth

We combine the output of our root emergence models and model 23 to
recreate root trajectories over time using `sim_plants.R` and
`sim_plants.stan`.
