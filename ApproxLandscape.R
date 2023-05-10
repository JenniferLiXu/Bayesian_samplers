#Define a mixture of two Gaussian distributions as our target T
target <- function(x, returnGrad = TRUE) {
  sample_1 = w * dnorm(x, mean = 2, sd = sigma1)
  sample_2 = (1 - w) * dnorm(x, mean = 2 - mu_target, sd = sigma2)
  distri = sample_1 + sample_2
  if (returnGrad) {
    grad = -(x-2)/sigma1 * sample_1 + -(x-2+mu_target)/sigma2  * sample_2
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
U <- function(q, returnGrad = TRUE) {
  distri = - log(target(q, returnGrad = FALSE))
  if(returnGrad){
    grad = - 1 / target(q, returnGrad = FALSE)* target(q)$grad
    return(list(distri = distri, grad = grad))
  }
  else return(distri) 
}

#initialize
set.seed(120)
n_iter <- 1000
#mu <- seq(10,18,length=4)
mu <- c(5)
sigma1 <- 0.5
sigma2 <- 1
w <- 0.3

# Run HMC for a short time
library(pracma)
source("HMC_MixGaussian.R")

mu_target <- mu
# Run HMC for the new target distribution T
samples <- hmc(U = U, epsilon = 0.1, L = 10, current_q = 0)
x <- samples$chain

par(mfrow = c(2,2))

#Print out the histogram of the samples
hist(x, breaks = 30, freq = FALSE, 
     main = substitute(paste("mixture gap mu=", a), 
                       list(a = mu_target) ))
U_x <- U(x)
plot(x, U_x$distri)

#Store the parameters of the target distribution
targetParams <- target(x)
#plot(x, targetParams$distri)

# Find q with the highest U
(highest_U_idx <- which.max(sapply(x, function(x) -U(x, returnGrad = FALSE))))
(q_with_highest_U <- x[highest_U_idx])

#Method "BFGS" is a quasi-Newton method (also known as a variable metric algorithm)
opt_result <- optim(q_with_highest_U, 
                    function(x) U(x, returnGrad = FALSE), 
                    function(x) U(x, returnGrad = TRUE)$grad, 
                    method ="BFGS",
                    hessian = TRUE)
(peak_U <- opt_result$value)
(peak_pos <- opt_result$par)
(hessian <- optimHess(peak_pos, function(x) U(x)$distri))

# Estimate the quadratic form using gradient of U near the peak
# the second derivative of U using finite difference method
(G1_mean <- peak_pos)
(G1_sd <- sqrt(1 / hessian))
g_total <- matrix(nrow = 1, ncol = 3)

g1 <- c(G1_mean, G1_sd, 1)

G <- function(x, g_info) {
  distri = g_info[3] * dnorm(x, mean = g_info[1], sd = g_info[2])
  grad = - (x - g_info[1])/g_info[2] * distri
  return(list(distri = distri, grad = grad))  
}

plot(x,G(x,g1)$distri)

# Define the new target distribution T' = T/G1
# and define the potential energy function U' and its gradient for T'
U_prime <- function(x ,returnGrad = TRUE) {
  U_x = U(x)
  G_x = G(x, g1)
  distri = U_x$distri + log(G_x$distri)
  if(returnGrad){
    grad = U_x$grad + G_x$grad/G_x$distri
    return(list(distri = distri, grad = grad))
  }
  else return(distri)
}

# Run HMC for the new target distribution T'
samples_prime <- hmc(U = U_prime, epsilon = 0.1, L = 10, current_q = 0)
x_prime <- samples_prime$chain
U_prime_x <- U_prime(x_prime)

hist(x_prime, breaks = 30, freq = FALSE, 
     main = substitute(paste("mixture gap mu=", a),
                       list(a = mu_target) ))
plot(x_prime, U_prime_x$distri)


# Find q with the highest U
(highest_U_idx_prime <- which.max(sapply(x_prime, function(x) - U_prime(x)$distri)))
(q_highest_U_prime <- x_prime[highest_U_idx_prime])

#Method "BFGS" is a quasi-Newton method (also known as a variable metric algorithm)
opt_result_prime <- optim(q_highest_U_prime, function(x) U_prime(x, returnGrad = FALSE), 
                    function(x) U_prime(x, returnGrad = TRUE)$grad, 
                    method ="BFGS",hessian = TRUE)
(peak_U_prime <- opt_result_prime$value)
(peak_pos_prime <- opt_result_prime$par)
(hessian_prime <- optimHess(peak_pos_prime, function(x) U_prime(x)$distri))

# Determine the position of the new peak and its quadratic form (second Gaussian G2)
(G2_mean <- peak_pos_prime)
(G2_sd <- sqrt(1 / hessian_prime))
 
g2 <- c(G2_mean, G2_sd, 1)

#Calculate the height of each peak
G1_height <- G(G1_mean,g1)$distri
G2_height <- G(G2_mean,g2)$distri

# Gaussian component information
# Gaussian component information
new_gaussian_components <- list(
  list(mean = G1_mean, sd = G1_sd, height = G1_height),
  list(mean = G2_mean, sd = G2_sd, height = G2_height)
)

# Calculate mixture weights
(total_height <- sum(sapply(new_gaussian_components, function(x) x$height)))
(mixture_weights <- sapply(new_gaussian_components, function(x) x$height / total_height))

# Define the mixture model M
M <- function(x) {
  distri <- 0
  grad = 0
  for (i in 1:length(gaussian_components)) {
    component <- gaussian_components[[i]]
    weight <- mixture_weights[i]
    info <- c(component$mean, component$sd, component$height)
    distri <- distri + G(x, info)$distri
    grad <- grad + G(x, info)$grad
  }
  return(list(distri = distri, grad = grad))
}

# Plot the recovered target distribution M
x_vals <- seq(-10, 10, 0.1)
recovered_vals <- sapply(x_vals, function(x) M(x)$distri)
plot(x_vals, recovered_vals, type = 'l', 
     main = "Recovered Target Distribution", 
     xlab = "x", ylab = "Density")

# Define the new target distribution T'' = T/M
# Define the potential energy function U'' and its gradient for T''
U_double_prime <- function(x ,returnGrad = TRUE) {
  U_x = U_prime(x)
  G_x = M(x)
  distri = U_x$distri + log(G_x$distri)
  if(returnGrad){
    grad = U_x$grad + G_x$grad/G_x$distri
    return(list(distri = distri, grad = grad))
  }
  else return(distri)
}


# Run HMC for the updated target distribution T''
n_iter <- 1000
samples_double_prime <- hmc(U = U_double_prime, epsilon = 0.1, L = 10, current_q = 0)
x_double_prime <- samples_double_prime$chain
hist(x_double_prime, breaks = 30, freq = FALSE, main = substitute(paste("mixture gap mu=", a), list(a = mu_target) ))
plot(U(x_double_prime)$distri)

# Evaluate the stopping criteria based on the mixture weights
#......

# Calculate the importance weights for the samples from M and (T/M)
importance_weights_M <- sapply(x, function(q) target(q)$distri / M(q))
importance_weights_T_over_M <- sapply(x_double_prime, function(q) M(q) / Target(q)$distri)
