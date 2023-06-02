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
  target_q = target(q)
  distri = - log(target_q$distri)
  if(returnGrad){
    grad = - 1 / target_q$distri* target_q$grad
    return(list(distri = distri, grad = grad))
  }
  else return(distri) 
}

#Find the value of the highest potential energy U and the mean and sd of the peak
peak <- function(x, U, M) {
  # Find q with the highest U
  highest_U_idx <- which.max(sapply(x, function(x) -U(x)$distri))
  q_with_highest_U <- x[highest_U_idx]
  
  #Method "BFGS" is a quasi-Newton method (also known as a variable metric algorithm)
  # Replace the target distribution T with a tempered distribution T' = T - M
  opt_result <- optim(q_with_highest_U, function(x) -target(x)$distri + M(x)$distri, 
                      function(x) -target(x)$grad + M(x)$grad, 
                      method ="BFGS",hessian = TRUE)
  peak_U <- opt_result$value
  peak_pos <- opt_result$par
  hessian <- optimHess(peak_pos, function(x) -target(x)$distri + M(x)$distri)
  
  # Estimate the quadratic form using gradient of U near the peak
  # the second derivative of U using finite difference method
  mean <- peak_pos
  sd <- sqrt(1 / hessian)
  
  return(list(peak_U = peak_U, mean = mean, sd = sd))
}

#Use g_info(=c(mean, sd, weight)) to generate a gaussian distribution function G
G <- function(x, g_info) {
  distri = g_info[3] * dnorm(x, mean = g_info[1], sd = g_info[2])
  grad = - (x - g_info[1])/g_info[2] * distri
  return(list(distri = distri, grad = grad))  
}

#Gaussian components, which will be put into a list
Gaussian_component <- function(peakU){
  # Estimate a Gaussian component
  G_component <- G(peakU$mean, c(peakU$mean, peakU$sd, 1))
  position = peakU$mean
  height <- G_component$distri
  
  return(list(mean = peakU$mean, sd = peakU$sd, height = height))
}

# Define the mixture model M
M <- function(x) {
  distri <- 0
  grad = 0
  for (i in 1:length(gaussians)) {
    component <- gaussians[[i]]
    weight <- mixture_weights[i]
    info <- c(component$mean, component$sd, weight)
    distri <- distri + G(x, info)$distri
    grad <- grad + G(x, info)$grad
  }
  return(list(distri = distri, grad = grad))
}

library(pracma)
source("HMC_MixGaussian.R")
par(mfrow = c(2,2))

max_iterations <- 3


#Initialize
set.seed(120)
n_iter <- 400
mu <- c(5)
sigma1 <- 0.5
sigma2 <- 1
w <- 0.3

mu_target <- mu

# Initialize the g_info and the list of Gaussian components
g_info <-c(0,0,0)
mixture_weights <- 1
components <-list(mean = 0, sd = 0, height = 0)
gaussians <- list(components)

# Run the loop
for (i in 1:2) {
  
  if (i == 1){
    UNew = U
  }
  else {
    # Adjusted UNew function by using U+log(M)
    UNew <- function(x, returnGrad = TRUE) {
      M_x = M(x)
      U_x = U(x)
      distri = U_x$distri + log(M_x$distri)
      if(returnGrad){
        grad = U_x$grad + M_x$grad/M_x$distri
        return(list(distri = distri, grad = grad))
      }
      else return(distri)
    }
  }

  # Run a short HMC chain to find a point with high potential energy U
  samples <- hmc(U = UNew, epsilon = 0.1, L = 10, current_q = 0)
  x <- samples$chain

  #Print out the histogram of the samples
  hist(x, breaks = 30, freq = FALSE, 
       main = substitute(paste("mixture gap mu=", a), list(a = mu_target)))
  plot(x, UNew(x)$distri)
  
  #find the peak of U and peak information
  peakU <- peak(x, UNew, M)
  g_info <- c(peakU$mean, peakU$sd, 1)
  print(g_info)
  
  components <- Gaussian_component(peakU)
  print(components)
  
  # Add the Gaussian component to the list
  gaussians[[i]] <- list(mean = components$mean, sd = components$sd, height = components$height)
  
  # Calculate mixture weights
  total_height <- sum(sapply(gaussians, function(x) x$height)) 
  mixture_weights <- sapply(gaussians, function(x) x$height / total_height)
  
  M = M
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

