#target T
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

#potential energy
U <- function(q, returnGrad = TRUE) {
  distri = - log(target(q, returnGrad = FALSE))
  if(returnGrad){
    grad = - 1 / target(q, returnGrad = FALSE)* target(q)$grad
    return(list(distri = distri, grad = grad))
  }
  else return(distri) 
}

# Initialize the list of Gaussian components
gaussians <- list()

peak <- function(x, U) {
  # Find q with the highest U
  highest_U_idx <- which.max(sapply(x, function(x) -U(x)$distri))
  q_with_highest_U <- x[highest_U_idx]
  
  #Method "BFGS" is a quasi-Newton method (also known as a variable metric algorithm)
  opt_result <- optim(q_with_highest_U, function(x) -target(x)$distri, 
                      function(x) -target(x)$grad, 
                      method ="BFGS",hessian = TRUE)
  peak_U <- opt_result$value
  peak_pos <- opt_result$par
  hessian <- optimHess(peak_pos, function(x) U(x)$distri)
  
  # Estimate the quadratic form using gradient of U near the peak
  # the second derivative of U using finite difference method
  mean <- peak_pos
  sd <- sqrt(1 / hessian)
  
  return(list(peak_U = peak_U, mean = mean, sd = sd))
}

G <- function(x, g_info) {
  distri = g_info[3] * dnorm(x, mean = g_info[1], sd = g_info[2])
  grad = - (x - g_info[1])/g_info[2] * distri
  return(list(distri = distri, grad = grad))  
}

Gaussian_component <- function(x, U, peakU, g_info){
  # Estimate a Gaussian component
  G_component <- G(x, g_info)
  height <- G(peakU$mean,g_info)$distri
  
  return(list(mean = peakU$mean, sd = peakU$sd, height = height))
}

# Define the mixture model M
M <- function(x) {
  distri <- 0
  grad = 0
  for (i in 1:length(gaussians)) {
    component <- gaussians[[i]]
    weight <- mixture_weights[i]
    info <- c(component$mean, component$sd, component$height)
    distri <- distri + G(x, info)$distri
    grad <- grad + G(x, info)$grad
  }
  return(list(distri = distri, grad = grad))
}


#initialize
set.seed(120)
n_iter <- 400
mu <- c(5)
sigma1 <- 0.5
sigma2 <- 1
w <- 0.3

mu_target <- mu

library(pracma)
source("HMC_MixGaussian.R")
par(mfrow = c(2,2))

max_iterations <- 3

# Run the loop
for (i in 1:max_iterations) {
  # Run a short HMC chain to find a point with high potential energy U
  samples <- hmc(U = U, epsilon = 0.1, L = 10, current_q = 0)
  x <- samples$chain
  
  #Print out the histogram of the samples
  hist(x, breaks = 30, freq = FALSE, 
       main = substitute(paste("mixture gap mu=", a), list(a = mu_target)))
  plot(x, U(x)$distri)
  
  peakU <- peak(x,U)
  g_info <- c(peakU$mean, peakU$sd, 1)
  
  components <- Gaussian_component(x, U, peakU, g_info)
  print(components)
  
  # Add the Gaussian component to the list
  gaussians[[i]] <- list(mean = components$mean, sd = components$sd, height = components$height)
  
  # Calculate mixture weights
  total_height <- sum(sapply(gaussians, function(x) x$height)) 
  mixture_weights <- sapply(gaussians, function(x) x$height / total_height)
  
  #Store the parameters of the target distribution
  targetParams <- target(x)
  
  G_x = G(x, g_info)
  
  # Replace the target distribution T with a tempered distribution T' = T / G
  # Define the potential energy function U' and its gradient for T'
  target <- function(x, returnGrad = TRUE){
    distri = targetParams$distri - G_x$distri
    grad = targetParams$grad - G_x$grad
    return(list(distri = distri, grad = grad))
  }
  
  U <- function(x ,returnGrad = TRUE) {
    U_x = U(x)
    distri = U_x$distri + log(G_x$distri)
    if(returnGrad){
      grad = U_x$grad + G_x$grad/G_x$distri
      return(list(distri = distri, grad = grad))
    }
    else return(distri)
  }
}


# Plot the recovered target distribution
x_vals <- seq(-5, 5, 0.1)
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

