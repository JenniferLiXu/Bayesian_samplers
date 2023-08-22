 # Repository Overview - Bayesian_samplers

This repository contains various files associated with our study comparing different sampling methods: MH-MCMC, HMC, and Iterative Gaussian Approximations for Hamiltonian Monte Carlo Sampling (IGA-HMC). Here's a brief rundown of what you can find in each file:

1. `MH_MCMC`: This file contains our implementation of the Metropolis-Hastings Markov Chain Monte Carlo (MH-MCMC) algorithm. The function `MH_MCMC(target, proposal, x_init, n_iter)` generates a Markov chain of samples from the target distribution using the proposal distribution. The function returns a list containing the sample chain and the maximum of the target function.

2. `HMC_MixGaussian`: This file contains our implementation of the Hamiltonian Monte Carlo (HMC) algorithm. The function `hmc(U, epsilon, L, current_q)` executes the HMC algorithm given a potential energy function, step size, number of steps, and initial position. The function returns a list containing the generated sample chain, acceptance rate, the position corresponding to the minimum potential energy, and the maximum of the target function.

3. `IGA_HMC_MixGaussian` : In this file we run the IGA-MCMC function to generate samples for fitting target distributions, which are a mixture of two Gaussians.

4. `result_plot`: In this file, we run the MH-MCMC and HMC functions to generate samples for fitting target distributions. These distributions are a mixture of two Gaussians with different distances between peaks.

5. `Approxlandscape`: This file introduces a novel method, Iterative Gaussian Approximations for Hamiltonian Monte Carlo Sampling (IGA-HMC). This method is used to generate samples for the mixture of Gaussians.

6. `landscape_loop`: In this file, we apply the IGA-HMC method with Importance Sampling on a real-world target distribution derived from a sum of Student’s t-distributions centered around the values of the `faithful$waiting` variable. 

7. `comparing(MH&HMC)`: In this final file, we compare the performances of the MH-MCMC and HMC algorithms when applied to the target distribution from the `landscape_loop` file.


