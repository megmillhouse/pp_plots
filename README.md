# pp_plots

A way to test that your MCMC sampler (or sampler of your choosing) is working properly is to check that for a large set of runs (say, 100) with parameters drawn from your prior, N% of the time the true (injected) value of the parameter will lie in the Nth credible interval. 

This code calculates the fraction of events where the parameter of interest is recovered within credible interval N as a function of N, and plots.
