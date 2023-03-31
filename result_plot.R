#Running and plotting of the results
#Define a mixture of two Gaussian distributions
target <- function(x, returnGrad = FALSE) {
  if(returnGrad){
    return(-(x-2)/sigma1 * w * dnorm(x, mean = 2, sd = sigma1)
           + -(x-2+mu_target)/sigma2 * (1 - w) * dnorm(x, mean = 2 - mu_target, sd = sigma2))
  }
  else return(w * dnorm(x, mean = 2, sd = sigma1) + (1 - w) * dnorm(x, mean = 2 - mu_target, sd = sigma2))
}

#potential energy
U <- function(q, returnGrad = FALSE) {
  if(returnGrad){
    return(1/target(q, returnGrad = FALSE) *target(q,returnGrad = TRUE) )
  }
  else return(-log(target(q, returnGrad = FALSE))) 
}

#proposal distribution function
proposal <- function(x) {
  rnorm(1, mean = x, sd = 0.5)
}

#initialize
set.seed(123)
n_iter <- 10000
mu <- seq(1,10,length=4)
sigma1 <- 1
sigma2 <- 2
w <- 0.6

for (j in 1:length(mu)){
  mu_target <- mu[j]
  
  #generate samples
  #using the Metropolis-Hastings algorithm
  #source("MH_MCMC.R")
  #x <- MH_MCMC(target, proposal, 0, n_iter) 

  #using Hamiltonian Monte Carlo
  source("HMC_MixGaussian.R")
  x <- hmc(U = U, epsilon = 0.25, L =5, current_q = 6)
  
  #target \int \sin(x) p(x) dx
  esti_sin <- mean(sin(x))
  
  #numerical integral
  num_sin <- integrate(function(x) sin(x) * target(x), -Inf, Inf)
  
  #difference between estimated value and numerical integral
  cat("difference:", abs(esti_sin - num_sin$value), "\n")
  #diff_est_num[j] <- abs(esti_sin - num_sin$value)
  # Print the acceptance rate
  #cat("Acceptance rate:", mean(accept), "\n")
  
  # Plot the samples and the target distribution
  hist(x, breaks = 30, freq = FALSE, main = substitute(paste("mixture gap mu=", a), list(a = mu_target) ))
  curve(target(x), add = TRUE, col = "red", lwd = 2)
  
  library(mixtools)
  #we are fitting a mixture of two Gaussian distributions
  fit <- normalmixEM(x, k = 2)
  # Plot the fitted mixture model
  curve(fit$lambda[1] * dnorm(x, fit$mu[1], sqrt(fit$sigma[1]^2)) +
          fit$lambda[2] * dnorm(x, fit$mu[2], sqrt(fit$sigma[2]^2)),
        from = min(x), to = max(x), add = TRUE, col = "blue")
  burnin <- 2000
  plot(x[-(1:burnin)], type = "l", xlab = "", main = "Chain values of x", )
}
