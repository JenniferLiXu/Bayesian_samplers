#initialize
set.seed(111)
n_iter <- 20000
mu <- seq(1,10,length=4)
#mu <- c(7)
sigma1 <- 1
sigma2 <- 2
w <- 0.7

diff_est_num <- numeric(length(mu))
par(mfrow = c(2,2))

#running and plotting of the results
#Define a mixture of two Gaussian distributions

target <- function(x, returnGrad = TRUE) {
  sample_1 = w * dnorm(x, mean = 2, sd = sigma1)
  sample_2 = (1 - w) * dnorm(x, mean = 2 - mu_target, sd = sigma2)
  distri = sample_1 + sample_2
  if (returnGrad) {
    grad = -(x-2)/sigma1^2 * sample_1 + -(x-2+mu_target)/sigma2^2  * sample_2
    return(list(distri = distri, grad = grad))
  } else {
    return(distri)
  }
}

#proposal distribution function
proposal <- function(x) {
  rnorm(1, mean = x, sd = 0.5)
}

#potential energy
#U <- function(q, returnGrad = FALSE) {
#  if(returnGrad){
#    return(1/target(q, returnGrad = FALSE) *target(q,returnGrad = TRUE) )
#  else return(-log(target(q, returnGrad = FALSE))) 
#  }
#}

U <- function(q, returnGrad = TRUE) {
  target_q = target(q)
  distri = - log(target_q$distri)
  if(returnGrad){
    grad = - 1 / target_q$distri* target_q$grad
    return(list(distri = distri, grad = grad))
  }
  else return(distri) 
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



for (j in 1:length(mu)){
  mu_target <- mu[j]
  
  library(R.utils)
  start_time <- Sys.time()
  source("MH_MCMC.R")
  #generate samples using the Metropolis-Hastings algorithm
  result <- MH_MCMC(target, proposal, 0, n_iter)
  
  #source("HMC_MixGaussian.R")
  #result <- hmc(U = U, epsilon = 0.25, L = 5, current_q = 0)
  
  end_time <- Sys.time()
  running_time <- end_time - start_time
  #Output the running time
  print(running_time)
  
  x <- result$chain
  
  #target \int \sin(x) p(x) dx
  esti_sin <- mean(sin(x))
  
  #numerical integral
  num_sin <- integrate(function(x) sin(x) * target(x)$distri, -Inf, Inf)
  
  #difference between estimated value and numerical integral
  cat("difference:", abs(esti_sin - num_sin$value), "\n")
  
  #diff_est_num[j] <- abs(esti_sin - num_sin$value)
  # Print the acceptance rate
  #cat("Acceptance rate:", mean(accept), "\n")
  
  # Plot the samples and the target distribution
  hist(x, breaks = 30, freq = FALSE, 
       main = substitute(paste("mixture gap mu=", a), list(a = mu_target)),
       ylim = c(0, result$f_max)
       )
  curve(target(x)$distri, add = TRUE, col = "red", lwd = 2)

  
  library(mixtools)
  #we are fitting a mixture of two Gaussian distributions
  fit <- normalmixEM(x, k = 2)
  # Plot the fitted mixture model
  curve(fit$lambda[1] * dnorm(x, fit$mu[1], sqrt(fit$sigma[1]^2)) +
          fit$lambda[2] * dnorm(x, fit$mu[2], sqrt(fit$sigma[2]^2)),
        from = min(x), to = max(x), add = TRUE, col = "blue")

  burnin <- 1000
  plot(x[-(1:burnin)], type = "l", xlab = "", main = "Chain values of x", )
  print(paste("The estimated weight w1:",chain_char(x[-(1:burnin)])$w_above))
}

