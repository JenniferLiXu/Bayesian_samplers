#Define a mixture of two Gaussian distributions
target <- function(x) {
  w * dnorm(x, mean = mu, sd = sigma1) + (1 - w) * dnorm(x, mean = 2 - mu_target, sd = sigma2)
}

#potential energy
U <- function(q) {
  -log(target(q))
}

grad_U <- function(q){
  1/target(q) *
  (w * (q-2)/sigma1*dnorm(q, mean = 2, sd = sigma1) 
  + (1 - w) * (q-2+mu_target)/sigma2* dnorm(q, mean = 2 - mu_target, sd = sigma2))
}


hmc <- function(U, grad_U, epsilon, L, current_q) {
  #input: 
  # U: the potential energy function, its gradient grad_U 
  # epsilon: the step size 
  # L: the number of steps
  # current_q: the current position 
  q <- current_q
  p <- matrix(rnorm(length(q),0,1), ncol = 1) 
  current_p <- p
  # Leapfrog method: 
  # First make a half step for momentum 
  p <- p - epsilon / 2 * grad_U(q)  
  for (i in 1:L) {
    # Make a full step for the position
    q <- q + epsilon * p 
    # Make a full step for the momentum, except at end of trajectory
    if (i != L) p <- p - epsilon * grad_U(q)
  }
  # Make a half step for momentum at the end
  p <- p - epsilon / 2 * grad_U(q)
  # Negate momentum at end of trajectory to make the proposal symmetric
  p <- -p
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U <- U(current_q)
  current_K <- sum(current_p^2) / 2
  proposed_U <- U(q)
  proposed_K <- sum(p^2) / 2
  
  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K)) {
    return(q)  #Accept
  } else {
    return(current_q)  #Reject
  }
  
}

set.seed(111)
n_iter <- 10000
#initialization
q_hmc <- numeric(n_iter)
mu <- seq(1,10,length=4)
q_init <- 6

sigma1 <- 1
sigma2 <- 2
w <- 0.6

par(mfrow = c(2,2))

for (j in 1:length(mu)){
  
  mu_target <- mu[j]
  
  for (i in 1:n_iter) {
    q_hmc[i] = hmc(U = U, grad_U = grad_U, epsilon = 0.25, L = 5, current = q_init)
    q_init = q_hmc[i]
  }
  burnin <- 1000
  x <- q_hmc[-(1:burnin)]
  hist(x, breaks = 50, freq = FALSE, main = substitute(paste("mixture gap mu=", a), list(a = mu_target) ))
  curve(target(x), add = TRUE, col = "red", lwd = 2)
  
  library(mixtools)
  #we are fitting a mixture of two Gaussian distributions
  fit <- normalmixEM(x, k = 2)
  # Plot the fitted mixture model
  curve(fit$lambda[1] * dnorm(x, fit$mu[1], sqrt(fit$sigma[1]^2)) +
          fit$lambda[2] * dnorm(x, fit$mu[2], sqrt(fit$sigma[2]^2)),
        from = min(x), to = max(x), add = TRUE, col = "blue")
  
  plot(x[-(1:burnin)], type = "l", xlab = "", main = "Chain values of q_hmc")
}