#running and plotting of the results
#Define a mixture of two Gaussian distributions
target <- function(x, returnGrad = FALSE) {
  if(returnGrad){
    return(-(x-2)/sigma1 * w * dnorm(x, mean = 2, sd = sigma1)
           + -(x-2+mu_target)/sigma2 * (1 - w) * dnorm(x, mean = 2 - mu_target, sd = sigma2))
  }
  else return(w * dnorm(x, mean = 2, sd = sigma1) + (1 - w) * dnorm(x, mean = 2 - mu_target, sd = sigma2))
}

#proposal distribution function
proposal <- function(x) {
  rnorm(1, mean = x, sd = 0.5)
}

#potential energy
U <- function(q, returnGrad = FALSE) {
  if(returnGrad){
    return(1/target(q, returnGrad = FALSE) *target(q,returnGrad = TRUE) )
  }
  else return(-log(target(q, returnGrad = FALSE))) 
}

#Counting the swaps and length of swaps to estimate the weights
chain_char <- function(chain){
  # Initialization
  average <- mean(chain)
  swaps <- 0
  swap_lengths <- c()
  current_length <- 0
  state <- integer(length(chain))
  state[1] <- ifelse(chain[1] >= average, "above", "below")
  
  for (i in 2:length(chain)) {
    state[i] <- ifelse(chain[i] >= average, "above", "below")
    if (state[i] != state[i-1]) {
      swaps <- swaps + 1
      swap_lengths[swaps] <- current_length
      current_length <- 1
    } else {
      current_length <- current_length + 1
    }
  }
  # Store the last swap length
  swap_lengths <- c(swap_lengths, current_length)
  
  # Calculate weights of being below or above average
  weights_below <- sum(state == "below") / length(chain)
  weights_above <- sum(state == "above") / length(chain)
  return(list(swaps = swaps, swap_lengths = swap_lengths, w_below = weights_below, 
              w_above = weights_above))
}

#initialize
set.seed(123)
n_iter <- 200
#mu <- seq(10,18,length=4)
mu <- c(7)
sigma1 <- 1
sigma2 <- 2
w <- 0.6

diff_est_num <- numeric(length(mu))
par(mfrow = c(2,2))

for (j in 1:length(mu)){
  
  mu_target <- mu[j]
  
  #source("~/Desktop/SemesterProject/Bayesian_samplers/MH_MCMC.R")
  #generate samples using the Metropolis-Hastings algorithm
  #x <- MH_MCMC(target, proposal, 0, n_iter)
  source("~/Desktop/SemesterProject/Bayesian_samplers/HMC_MixGaussian.R")
  result <- hmc(U = U, epsilon = 0.5, L = 5, current_q = 0)
  x <- result$chain
  #target \int \sin(x) p(x) dx
  esti_sin <- mean(sin(x))
  
  #numerical integral
  num_sin <- integrate(function(x) sin(x) * target(x), -Inf, Inf)
  
  #difference between estimated value and numerical integral
  #cat("difference:", abs(esti_sin - num_sin$value), "\n")
  
  #diff_est_num[j] <- abs(esti_sin - num_sin$value)
  # Print the acceptance rate
  #cat("Acceptance rate:", mean(accept), "\n")
  
  # Plot the samples and the target distribution
  hist(x, breaks = 30, freq = FALSE, main = substitute(paste("mixture gap mu=", a), list(a = mu_target) ))
  curve(target(x), add = TRUE, col = "red", lwd = 2)
  
  library(mixtools)
  #we are fitting a mixture of two Gaussian distributions
  #fit <- normalmixEM(x, k = 2)
  # Plot the fitted mixture model
  curve(fit$lambda[1] * dnorm(x, fit$mu[1], sqrt(fit$sigma[1]^2)) +
          fit$lambda[2] * dnorm(x, fit$mu[2], sqrt(fit$sigma[2]^2)),
        from = min(x), to = max(x), add = TRUE, col = "blue")

  burnin <- 100
  plot(x[-(1:burnin)], type = "l", xlab = "", main = "Chain values of x", )
  print(paste("The estimated weight w1:",chain_char(x[-(1:burnin)])$w_above))
}

