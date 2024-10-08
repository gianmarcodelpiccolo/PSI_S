Stochastic trans-dimensional code for seismic imaging.

Developed by G. Del Piccolo, PhD student at the Department of Geosciences, University of Padova.

Published reference for the method: G. Del Piccolo, B.P. VanderBeek, M. Faccenda, A. Morelli, J.S. Byrnes; Imaging Upper‐Mantle Anisotropy with Transdimensional Bayesian Monte Carlo Sampling. Bulletin of the Seismological Society of America 2024; 114 (3): 1214–1226. doi: https://doi.org/10.1785/0120230233

Build the inversion file running:

    julia buildIP.jl

To run the solver:

    julia run_rjMCMC.jl $(run_name) $(number_of_chains) $(batch_number)

Merge the chains:

    julia merge_chains.jl $(run_name) $(number_of_chains) $(discarded_models_per_chain)

Some plotting scripts are available in the tutorial folders.


