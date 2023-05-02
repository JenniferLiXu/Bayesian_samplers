#Define a mixture of two Gaussian distributions as our target T
target <- function(x) {
  sample_1 = w * dnorm(x, mean = 2, sd = sigma1)
  sample_2 = (1 - w) * dnorm(x, mean = 2 - mu_target, sd = sigma2)
  distri = sample_1 + sample_2
  grad = -(x-2)/sigma1 * sample_1 + -(x-2+mu_target)/sigma2  * sample_2
  return(list(distri = distri, grad = grad))
}

#proposal distribution function
proposal <- function(x) {
  rnorm(1, mean = x, sd = 0.5)
}

#potential energy
U <- function(q) {
  #Return U
  target_q = target(q)
  distri = -log(target_q$distri)
  #Return dU
  grad = 1/target_q$distri * target_q$grad
  #grad = (-log(target(q + 1e-5)$distri) + log(target(q - 1e-5)$distri))/ (2 * 1e-5)
  return(list(distri = distri, grad = grad))
}

#initialize
set.seed(23)
n_iter <- 1000
#mu <- seq(10,18,length=4)
mu <- c(5)
sigma1 <- 0.5
sigma2 <- 1
w <- 0.5

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
plot(x, U(x)$distr)
plot(x, target(x)$distri)

# Find q with the highest U
(highest_U_idx <- which.max(sapply(x, function(x) -U(x)$distri)))
(q_with_highest_U <- x[highest_U_idx])

# Perform numerical optimization to find the peak of U
# opt_result <- steep_descent(q_with_highest_U, function(x) U(x)$distri, function(x) U(x)$grad)
# (peak_U <- opt_result$fmin)
# (peak_pos <- opt_result$xmin)

#Method "BFGS" is a quasi-Newton method (also known as a variable metric algorithm)
opt_result <- optim(q_with_highest_U, function(x) U(x)$distri, function(x) U(x)$grad, 
                    method ="BFGS",hessian = TRUE)
(peak_U <- opt_result$value)
(peak_pos <- opt_result$par)
(hessian <- optimHess(peak_pos, function(x) U(x)$distri))

#opt_result <- optimize(function(x) U(x)$distri, interval = range(x), maximum = FALSE)
#peak_U <- opt_result$objective
#peak_pos <- opt_result$minimum

# Estimate the quadratic form using gradient of U near the peak
# the second derivative of U using finite difference method
#(hessian <- (U(peak_pos + 1e-5)$grad  - U(peak_pos - 1e-5)$grad ) / (2 * 1e-5))
(G1_mean <- peak_pos)
(G1_sd <- sqrt(1 / hessian))

# Print the results
cat("Estimated first Gaussian component G1:\n")
cat(paste("Mean:", peak_pos, "\n"))


G1 <- function(x) {
  dnorm(x, mean = G1_mean, sd = G1_sd)
}

plot(x,G1(x))

# Define the new target distribution T' = T/G1
new_target <- function(x) {
  target_x <- target(x)
  distri = target_x$distri / G1(x) 
  #grad = (U(peak_pos + 1e-5)$grad  - U(peak_pos - 1e-5)$grad ) / (2 * 1e-5)
  grad = target_x$grad / G1(x)
  return(list(distri = distri, grad = grad))
}


# Define the potential energy function U' and its gradient for T'
U_prime <- function(x) {
  new_target_x = new_target(x)
  target_x = target(x)$distri
  #distri = -log(new_target_x$distri)
  distri = - log(target_x) + log(G1(x))
  #grad = (-log(new_target(x + 1e-5)$distri) + log(new_target(x + 1e-5)$distri) ) / (2 * 1e-5)
  grad = 1/new_target_x$distri * new_target_x$grad
  return(list(distri = distri, grad = grad))
}

plot(x, U_prime(x)$distr)

# Run HMC for the new target distribution T'
#samples_prime <- hmc(U = U_prime, epsilon = 0.1, L = 10, current_q = 0)
samples_prime <- hmc(U = U_prime, epsilon = 0.1, L = 10, current_q = 0)
x_prime <- samples_prime$chain

hist(x_prime, breaks = 30, freq = FALSE, main = substitute(paste("mixture gap mu=", a), list(a = mu_target) ))
plot(x_prime, U(x_prime)$distri - log(G1(x_prime)))
plot(x_prime,new_target(x)$distri)


# Find q with the highest U
(highest_U_idx_prime <- which.max(sapply(x_prime, function(x) - U_prime(x)$distri)))
(q_highest_U_prime <- x_prime[highest_U_idx_prime])

#Method "BFGS" is a quasi-Newton method (also known as a variable metric algorithm)
opt_result_prime <- optim(q_highest_U_prime, function(x) U_prime(x)$distri, 
                    function(x) U_prime(x)$grad, 
                    method ="BFGS",hessian = TRUE)
(peak_U_prime <- opt_result_prime$value)
(peak_pos_prime <- opt_result_prime$par)
(hessian_prime <- optimHess(peak_pos, function(x) U_prime(x)$distri))

# Perform numerical optimization (gradient decent) to find the peak of U
# opt_result_prime <- steep_descent(q_highest_U_prime, function(x) U_prime(x)$distri,
#                                   function(x) U_prime(x)$grad)
# (peak_U_prime <- -opt_result_prime$fmin)
# (peak_pos_prime <- opt_result_prime$xmin)

# Estimate the quadratic form using gradient of U' near the new peak
# the second derivative of U using finite difference method
# hessian_prime <- (U_prime(peak_pos + 1e-5)$grad  - U_prime(peak_pos - 1e-5)$grad ) / (2 * 1e-5)

# Determine the position of the new peak and its quadratic form (second Gaussian G2)G2_mean <- peak_pos_prime
G2_mean <- peak_pos_prime
G2_sd <- sqrt(1 / hessian_prime)

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

plot(x,M(x))

# Define the new target distribution T'' = T/M
target_double_prime <- function(x) {
  distri = target(x)$distri / M(x)
  grad = target(x)$grad / M(x)
  return(list(distri = distri, grad = grad))
}

plot(x, target_double_prime(x)$distri)

# Define the potential energy function U'' and its gradient for T''
U_double_prime <- function(x) {
  target_double_prime_x = target_double_prime(x)
  distri = -log(target_double_prime_x$distri)
  grad =  1/target_double_prime_x$distri * target_double_prime_x$grad
  return(list(distri = distri, grad = grad))
}

# Run HMC for the updated target distribution T''
set.seed(42)
n_samples <- 1000
samples_double_prime <- hmc(U = U_double_prime, epsilon = 0.1, L = 10, current_q = 0)
x_double_prime <- samples_double_prime$chain
hist(x_double_prime, breaks = 30, freq = FALSE, main = substitute(paste("mixture gap mu=", a), list(a = mu_target) ))
plot(U(x_double_prime)$distri)
