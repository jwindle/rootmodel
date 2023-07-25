# Overview

To reproduce the plots in the "Identifiability" section of the paper
walk through `sim_plants_v1.R`.  That file has two broad parts.
First, it will create a factorial combination of parameter values and
then simulate data for each parameter combination.  Second, it will
make a series of plots, some of which are used in the paper to provide
the numerical justification for a lack of identifiable parameters.

If you want to see the impact of that on the posterior, then walk
through `fit_plants_v1.R`.  That file will try to fit the model to one
of the simulated data sets.  The sampler is slow and it shows how you
do not concentrate on a point, but on a surface.


