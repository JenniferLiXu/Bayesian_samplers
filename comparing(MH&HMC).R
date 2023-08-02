set.seed(123)
n_iter <- 400

par(mfrow = c(1,1))
#Example (unnormalised) target:
likli <- function(x) {
  # for df = 1, we have f(z) = log (1 / (π * (1 + z²)) =  -log(π) - log(1 + z²)
  # sum(f(z))/500 
  
  #for df = 1, we have
  sum(dt(x - faithful$waiting, 1, log = TRUE))/500
}

mus <- 400:1000/10
lls <- sapply(mus, likli)

plot(mus, exp(lls))

unnormalize_target <- function(x, returnGrad = TRUE){
  lls = sapply(x, likli)
  distri = exp(lls)
  if(returnGrad) {
    grad = sapply(x, function(x_i) {
      z = x_i - faithful$waiting
      distri_i = exp(likli(x_i))
      #for df = 1
      distri_i * sum((-2* z/(1 + z^2)))/500
      
      #for df = 3
      #distri_i * (-4 /3) * sum( z/(1 + (z^2)/3))/500
    }) 
    return(list(distri = distri, grad = grad))
  } else{
    return(distri)
  }
}
# Calculate the total area under the curve by integrating over the entire range
total_area <- integrate(Vectorize(function(x) unnormalize_target(x, returnGrad = FALSE)), lower = 40, upper = 100)$value

target <- function(x, returnGrad = TRUE){
  unnormalize_target_x = unnormalize_target(x)
  # normalize the distribution
  normalized_distri = unnormalize_target_x$distri / total_area
  if(returnGrad) {
    # normalize the gradient
    normalized_distri = normalized_distri
    normalized_grad = unnormalize_target_x$grad / total_area
    
    return(list(distri = normalized_distri, grad = normalized_grad))
  } else{
    # normalize only the distribution
    return(normalized_distri)
  }
}

#proposal distribution function
proposal <- function(x) {
  rnorm(1, mean = x, sd = 0.5)
}


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

library(R.utils)


start_time <- Sys.time()
#withTimeout({
  
  source("MH_MCMC.R")
  #generate samples using the Metropolis-Hastings algorithm
  result <- MH_MCMC(target, proposal, 50, n_iter)
  
  #source("HMC_MixGaussian.R")
  #result <- hmc(U = U, epsilon = 0.25, L = 5, current_q = 50)

#}, timeout = 5, onTimeout = "silent")  # Time limit of 5s

end_time <- Sys.time()
  
running_time <- end_time - start_time

print(running_time)
x <- result$chain



# Plot the samples and the target distribution
hist(x, breaks = 30, freq = FALSE, 
     main = "Histogram of samples(MH-MCMC)",
     #main = "Histogram of samples(HMC)",
     ylim = c(0, 0.035)
)
x_vals <- seq(20, 120, 0.1)
lines(x_vals, target(x_vals)$distri, type = 'l', col = "black")

library(mixtools)
#we are fitting a mixture of two Gaussian distributions
fit <- normalmixEM(x, k = 2)
# Plot the fitted mixture model
curve(fit$lambda[1] * dnorm(x, fit$mu[1], sqrt(fit$sigma[1]^2)) +
        fit$lambda[2] * dnorm(x, fit$mu[2], sqrt(fit$sigma[2]^2)),
      from = min(x), to = max(x), add = TRUE, col = "blue")

# Add legend
legend("topright", # place it at the top right corner
       legend = c("Target Distribution", "Recovered Distribution"), 
       col = c("black", "blue"), # colors should be in the same order
       lty = 1, # line types
       cex = 0.5) # control size of legend

print(fit$mu)

