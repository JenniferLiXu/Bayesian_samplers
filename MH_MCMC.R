#initialize
set.seed(123)
n <- 10000
mu <- seq(1,10,length=4)
sigma1 <- 1
sigma2 <- 2
w <- 0.6
sd_proposal <- 0.5

#proposal distribution function
proposal <- function(x) {
  rnorm(1, mean = x, sd = sd_proposal)
}

accept <- numeric(n)
diff_est_num <- numeric(length(mu))

MH_MCMC <- function(target, proposal, x_init, n_iter) {
  #target: the target distribution
  #proposal: the proposal distribution function
  #x_init: starting value of the chain
  #n_iter: number of iterations
  chain <- numeric(n_iter)
  chain[1] <- x_init
  f_value <- target(x_init)
  
  for(i in 2:n_iter){
    # Generate candidate sample
    y <- proposal(chain[i - 1])
    f_proposed <- target(y)
    # Calculate acceptance ratio
    alpha <- min(1, f_proposed / f_value)
    
    # Decide whether to accept or reject the candidate sample
    if (runif(1) < alpha) {
      chain[i] <- y
      accept[i] <- 1
      f_value <- f_proposed
    }else {
      chain[i] <- chain[i - 1]
      accept[i] <- 0
    }
  }
  return(chain)
}

par(mfrow = c(2,2))
  
for (j in 1:length(mu)){
  
  mu_target <- mu[j]
  #Define a mixture of two Gaussian distributions
  target <- function(x) {
    w * dnorm(x, mean = 2, sd = sigma1) + (1 - w) * dnorm(x, mean = 2 - mu_target, sd = sigma2)
  }

  #generate samples using the Metropolis-Hastings algorithm
  x <- MH_MCMC(target, proposal, 0, n)
  
  #target \int \sin(x) p(x) dx
  esti_sin <- mean(sin(x))
  
  #numerical integral
  num_sin <- integrate(function(x) sin(x) * target(x), -Inf, Inf)
  
  #difference between estimated value and numerical integral
  cat("difference:", abs(esti_sin - num_sin$value), "\n")
  #diff_est_num[j] <- abs(esti_sin - num_sin$value)
  # Print the acceptance rate
  #cat("Acceptance rate:", mean(accept), "\n")
  
  # Plot the samples and the target distribution
  hist(x, breaks = 30, freq = FALSE, main = substitute(paste("mixture gap mu=", a), list(a = mu_target) ))
  curve(target(x), add = TRUE, col = "red", lwd = 2)
  
  library(mixtools)
  #we are fitting a mixture of two Gaussian distributions
  fit <- normalmixEM(x, k = 2)
  # Plot the fitted mixture model
  curve(fit$lambda[1] * dnorm(x, fit$mu[1], sqrt(fit$sigma[1]^2)) +
          fit$lambda[2] * dnorm(x, fit$mu[2], sqrt(fit$sigma[2]^2)),
        from = min(x), to = max(x), add = TRUE, col = "blue")
  burnin <- 2000
  plot(x[-(1:burnin)], type = "l", xlab = "", main = "Chain values of x", )
}
