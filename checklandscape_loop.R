#Initialize
set.seed(123)
n_iter <- 400
start_point <- 50

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

#Find the value of the highest T_prime and estimate a Gaussian component
estimate_gaussian <- function(peak_pos, M , T) {
  #Method "BFGS" is a quasi-Newton method (also known as a variable metric algorithm)
  #Replace the target distribution T with a tempered distribution T' = T - M
  opt_result <- optim(peak_pos, function(x) -T(x)$distri, function(x) -T(x)$grad, 
                      method ="BFGS",hessian = TRUE)
  peak_U <- opt_result$value
  peak_pos_T <- opt_result$par
  #hessian <- optimHess(peak_pos, function(x) -U(x)$distri, function(x) -U(x)$grad) #Compute hessian using U
  hessian <- optimHess(peak_pos_T, function(x) target(x)$distri, function(x) target(x)$grad)
  
  print(hessian)
  
  # Estimate the quadratic form using gradient of U near the peak
  # the second derivative of U using finite difference method
  mean <- peak_pos_T
  sd <- (-target(peak_pos_T)$distri/(hessian))^(1/2)
  
  return(list(peak_U = peak_U, mean = mean, sd = sd, hessian = hessian))
}

#Use g_info(=c(mean, sd, weight)) to generate a gaussian distribution function
gaussian_pdf <- function(x, g_info) {
  distri = g_info[3] * dnorm(x, mean = g_info[1], sd = g_info[2])
  grad = - (x - g_info[1])/g_info[2] * distri
  #distri = g_info[3] * exp(-(x - g_info[1])^2 / (2 * g_info[2]^2))
  #grad = - (x - g_info[1])/g_info[2] * distri
  return(list(distri = distri, grad = grad))  
}  

# Define the mixture model M
M <- function(x) {
  distri <- 0
  grad = 0
  for (i in 1:length(gaussians)) {
    component <- gaussians[[i]]
    weight_i <- mixture_weights[i]
    info <- c(component$mean, component$sd, weight_i)
    g_new <- gaussian_pdf(x, info)
    distri <- distri + g_new$distri
    grad <- grad + g_new$grad
  }
  return(list(distri = distri, grad = grad))
}

#Update target distribution 
T_prime <- function(x, returnGrad = TRUE){
  target_x = target(x)
  distri = target_x$distri - M(x)$distri
  #distri_log = log(distri)
  if (returnGrad) {
    grad = target_x$grad - M(x)$grad
    #grad_log = 1/distri * grad
    #return(list(distri = distri, grad = grad, distri_log = distri_log, grad_log = grad_log))
    return(list(distri = distri, grad = grad))
  } else {
    return(distri)
  }
}


library(pracma)
source("HMC_MixGaussian.R")

max_iterations <- 15

# Initialize the g_info and the list of Gaussian components
g_info <-c(0,0,0)
mixture_weights <- 1
components <-list(mean = 0, sd = 0)
gaussians <- list(components)
total_T <- c()


# Define a threshold for the weight below which the loop should stop
weight_threshold <- 0.02
weight_diff_threshold <- 0.003

# Small regularization term in the denominator to avoid extremely large values
reg_epsilon <- 0

library(R.utils)
start_time <- Sys.time()
# Run the loop
#withTimeout({
for (i in 1:max_iterations) {
  if (i == 1){
    Unew = U
    T = target
  }else {
    # Adjusted UNew function by using U+log(M)
    U_prime <- function(x, returnGrad = TRUE) {
      M_x = M(x)
      target_x = target(x)
      U_x = U(x)
      distri = -log(target_x$distri/(M_x$distri+ reg_epsilon))
      # distri = -log(target_x$distri - 0.02*M_x$distri)
      # distri = U_x$distri + log(M_x$distri)
      if(returnGrad){
        grad = U_x$grad + M_x$grad/(M_x$distri+reg_epsilon)
        # grad = -(target_x$grad-0.02*M_x$grad)/(target_x$distri - 0.02*M_x$distri)
        return(list(distri = distri, grad = grad))
      }
      else return(distri)
    }
    Unew = U_prime
    # Plot the recovered U
    x_vals <- seq(40, 100, 0.1)
    recovered_vals <- sapply(x_vals, function(x) Unew(x)$distri)
    plot(x_vals, recovered_vals, type = 'l', 
         main = "Recovered U", 
         xlab = "x", ylab = "Density")
  }
  #Run a short HMC chain to find a point with high potential energy U
  samples <- hmc(U = Unew, epsilon = 0.1, L = 10, current_q = start_point)
  x <- samples$chain
  # Find the peak of U directly from the HMC chain
  (min_q <- samples$min_q)
  
  #Print out the histogram of the samples
  #hist(x, breaks = 30, freq = FALSE, 
  #     main = substitute(paste("mixture gap mu=", a), list(a = mu_target)))
  hist(x, breaks = 30, freq = FALSE, main = "Histogram of samples")
  Unew_x = sapply(x, function(x) Unew(x)$distri)
  plot(x, Unew_x)
  
  #find the estimate_gaussian of U and peak information
  peakU <- estimate_gaussian(peak_pos = min_q, M, T)
  g_info <- c(peakU$mean, peakU$sd, 1)
  print(g_info)
  
  # If peakU$sd is NaN, stop the loop
  if (is.nan(peakU$sd)) {
    print("peakU$sd is NaN, stopping the loop.")
    break
  }
  
  # Calculate total likelihood for the new component
  new_T <- T(peakU$mean)$distri* peakU$sd
  
  # Compute new component's weight
  new_weight <- new_T / (sum(total_T) + new_T)
  
  # Check if the weight of the current component is below the threshold
  if (new_weight < weight_threshold || T_prime(peakU$mean)$distri/ (sum(total_T) + new_T) < weight_diff_threshold) {
    print(paste("Weight of Gaussian component", i, "is below the threshold. Stopping the loop."))
    break
  } else {
    # Add the new likelihood to the total_likelihoods list
    total_T <- c(total_T, new_T)
    
    # Recalculate weights: old weights stay the same relative to each other, new weight is calculated relative to total
    mixture_weights <- total_T/sum(total_T)
    
    # Add the Gaussian component to the list
    gaussians[[i]] <- list(mean = peakU$mean, sd = peakU$sd)
  } 
  
  M = M
  #plot(x, log(M(x)$distri))
  T = T_prime
  
}
#}, timeout = 5)  # Time limit of 5 second
end_time <- Sys.time()

running_time <- end_time - start_time

print(running_time)


# Plot the recovered target distribution
x_vals <- seq(20, 120, 0.1)
recovered_vals <- sapply(x_vals, function(x) gaussian_pdf(x, g_info)$distri)
recovered_vals <- sapply(x_vals, function(x) M(x)$distri)
plot(x_vals, target(x_vals)$distri, type = 'l', 
     main = "Recovered Target Distribution", 
     ylim = c(0, 0.06),
     xlab = "x", ylab = "Density", col = "black")
lines(x_vals, recovered_vals, col = "blue")

# Add legend
legend("topright", # place it at the top right corner
       legend = c("Target Distribution", "Recovered Distribution"), 
       col = c("black", "blue"), # colors should be in the same order
       lty = 1, # line types
       cex = 0.5) # control size of legend


##Importance sampling

# Function to sample from the mixture
sample_from_mixture <- function(n) {
  # Determine how many samples to draw from each Gaussian
  n_per_gaussian <- round(n * mixture_weights)
  
  # Ensure we get 'n' samples
  n_per_gaussian[1] <- n - sum(n_per_gaussian[-1])
  
  samples <- numeric(n)
  index <- 1
  
  for (i in 1:length(gaussians)) {
    component <- gaussians[[i]]
    samples[index:(index + n_per_gaussian[i] - 1)] <- 
      rnorm(n_per_gaussian[i], mean=component$mean, sd=component$sd)
    
    index <- index + n_per_gaussian[i]
  }
  
  return(samples)
}

# Draw 1000 samples
is_samples <- sample_from_mixture(2000)

# Compute importance weights
is_weights <- numeric(2000)

for(i in 1:2000) {
  is_weights[i] <- target(is_samples[i])$distri / M(is_samples[i])$distri
}

# Normalize the weights
is_weights <- is_weights / sum(is_weights)

# Use the weighted samples for inference
# target \int \sin(x) p(x) dx
esti_sin <- sum(is_weights * sin(is_samples))
# numerical integral
num_sin <- integrate(function(x) sin(x) * target(x)$distri, 40, 100)
# difference between estimated value and numerical integral
cat("difference:", abs(esti_sin - num_sin$value), "\n")

# Plot the samples
x <- is_samples
hist(is_samples, breaks = 50, freq = FALSE, main = "Histogram of samples")
curve(target(x)$distri, add = TRUE, col = "red", lwd = 2)

# Add legend
legend("topright", # place it at the top right corner
       legend = c("Target Distribution", "Recovered Distribution"), 
       col = c("red", "blue"), # colors should be in the same order
       lty = 1, # line types
       cex = 0.35) # control size of legend

# resampled_indices <- sample(1:length(is_samples), length(is_samples), 
#                             replace=TRUE, prob=is_weights)
# resampled_samples <- is_samples[resampled_indices]
# x <- resampled_samples
# hist(is_samples, breaks = 50, freq = FALSE, main = "Histogram of resampled samples")
