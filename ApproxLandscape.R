#Define a mixture of two Gaussian distributions as our target T
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

#initialize
set.seed(123)
n_iter <- 1000
#mu <- seq(10,18,length=4)
mu <- c(7)
sigma1 <- 1
sigma2 <- 2
w <- 0.6

# Run HMC for a short time
library(pracma)
#q0 <- rnorm(1)
source("~/Desktop/SemesterProject/Bayesian_samplers/HMC_MixGaussian.R")

mu_target <- mu
samples <- hmc(U = U, epsilon = 0.1, L = 10, current_q = 0)
x <- samples$chain

par(mfrow = c(2,2))
#Print out the histogram of the samples
hist(x, breaks = 30, freq = FALSE, main = substitute(paste("mixture gap mu=", a), list(a = mu_target) ))
plot(U(x))

# Find q with the highest U
highest_U_idx <- which.max(sapply(x, U))
q_with_highest_U <- x[highest_U_idx]

# Perform numerical optimization (gradient decent) to find the peak of U
opt_result <- steep_descent(q_with_highest_U, function(x) -U(x), function(x) -U(x, returnGrad = TRUE))
peak_U <- -opt_result$fmin
peak_pos <- opt_result$xmin

# Find the peak of U (minimum of -U)
#opt_result <- optimize(function(x) -U(x), interval = range(x))
#peak_pos <- opt_result$minimum
#peak_U <- -opt_result$objective

# Estimate the quadratic form using gradient of U near the peak
# the second derivative of U using finite difference method
hessian <- (U(peak_pos + 1e-5, returnGrad = TRUE)  - U(peak_pos - 1e-5, returnGrad = TRUE) ) / (2 * 1e-5)
G1_mean <- peak_pos
G1_sd <- sqrt(-1 / hessian)

# Print the results
cat("Estimated first Gaussian component G1:\n")
cat(paste("Mean:", peak_pos, "\n"))


G1 <- function(x) {
  dnorm(x, mean = G1_mean, sd = G1_sd)
}

# Define the new target distribution T' = T/G1
new_target <- function(x, returnGrad = FALSE) {
  if(returnGrad){
    return(target(x, returnGrad = TRUE) / G1(x) )
  }
  else return (target(x) / G1(x) )
}

# Define the potential energy function U' and its gradient for T'
'''
U_prime <- function(x, returnGrad = FALSE) {
  if(returnGrad){
    return(1/new_target(x) * new_target(x, returnGrad = TRUE) )
  }
  else return(-log(new_target(x)))
}
'''

U_prime <- function(x, returnGrad = FALSE) {
  if(returnGrad){
    derivative <- (log(new_target(x - 1e-5)) -log(new_target(x + 1e-5)) ) / (2 * 1e-5)
    return(derivative)       
  }
  else return(-log(new_target(x)))
}


# Run HMC for the new target distribution T'
#samples_prime <- hmc(U = U_prime, epsilon = 0.1, L = 10, current_q = 0)
samples_prime <- hmc(U = U_prime, epsilon = 0.1, L = 10, current_q = 0)
x_prime <- samples_prime$chain
hist(x_prime, breaks = 30, freq = FALSE, main = substitute(paste("mixture gap mu=", a), list(a = mu_target) ))
plot(x_prime)
# Find q with the highest U
highest_U_idx_prime <- which.max(sapply(x_prime, U_prime))
q_with_highest_U_prime <- x[highest_U_idx_prime]

# Perform numerical optimization (gradient decent) to find the peak of U
opt_result_prime <- steep_descent(q_with_highest_U_prime, function(x) -U_prime(x), function(x) -U_prime(x, returnGrad = TRUE))
(peak_U_prime <- -opt_result_prime$fmin)
(peak_pos_prime <- opt_result_prime$xmin)

# Estimate the quadratic form using gradient of U' near the new peak
# the second derivative of U using finite difference method
hessian_prime <- (U_prime(peak_pos + 1e-5, returnGrad = TRUE)  - U_prime(peak_pos - 1e-5, returnGrad = TRUE) ) / (2 * 1e-5)

# Determine the position of the new peak and its quadratic form (second Gaussian G2)G2_mean <- peak_pos_prime
G2_mean <- peak_pos_prime
G2_sd <- sqrt(-1 / hessian_prime)

cat("Estimated second Gaussian component (G2):\n")
cat(paste("Mean:", G2_mean, "\n"))
cat(paste("Standard deviation:", G2_sd, "\n"))
 

# Calculate the mixture weights
nu_1 <- dnorm(x, G1_mean, G1_sd) / (dnorm(x, G1_mean, G1_sd) + dnorm(x, G2_mean, G2_sd))
nu_2 <- 1 - nu_1

# Define the mixture model M
M <- function(x) {
  nu_1 * dnorm(x, G1_mean, G1_sd) + nu_2 * dnorm(x, G2_mean, G2_sd)
}

# Define the new target distribution T'' = T/M
target_double_prime <- function(x, returnGrad = FALSE) {
  if(returnGrad){
    target(x, returnGrad = TRUE) / M(x)
  }
  else target(x) / M(x)
}

# Define the potential energy function U'' and its gradient for T''
U_double_prime <- function(x, returnGrad = FALSE) {
  if(returnGrad){
    return(1/target_double_prime(x) * target_double_prime(x, returnGrad = TRUE) )
  }
  else return(-log(target_double_prime(x)))
}

# Run HMC for the updated target distribution T''
set.seed(42)
init_q <- 0
n_samples <- 100
samples_double_prime <- hmc(U = U_double_prime, epsilon = 0.1, L = 10, current_q = init_q)
x_double_prime <- samples_double_prime$chain
hist(x_double_prime, breaks = 30, freq = FALSE, main = substitute(paste("mixture gap mu=", a), list(a = mu_target) ))
plot(U(x_double_prime))
