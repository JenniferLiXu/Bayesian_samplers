#Define a mixture of two Gaussian distributions
f <- function(x, mu1, mu2, sigma1, sigma2, w) {
  w * dnorm(x, mean = mu1, sd = sigma1) + (1 - w) * dnorm(x, mean = mu2, sd = sigma2)
}

#proposal distribution function
proposal <- function(x, sd) {
  rnorm(1, mean = x, sd = sd)
}

#initialize
set.seed(123)
n <- 1000
x <- numeric(n)
accept <- numeric(n)

mu1 <- 1
mu2 <- -1
sigma1 <- 1
sigma2 <- 2
w <- 0.5
sd = 0.5

#generate samples using the Metropolis-Hastings algorithm
for (i in 2:n) {
  # Generate candidate sample
  y <- proposal(x[i-1], sd)
  
  # Calculate acceptance ratio
  alpha <- min(1, f(y, mu1, mu2, sigma1, sigma2, w) / f(x[i-1], mu1, mu2, sigma1, sigma2, w) * dnorm(x[i-1], mean = y, sd = sd) / dnorm(y, mean = x[i-1], sd = sd))
  
  # Decide whether to accept or reject the candidate sample
  if (runif(1) < alpha) {
    x[i] <- y
    accept[i] <- 1
  }else {
    x[i] <- x[i-1]
    accept[i] <- 0
  }
  

}

# Print the acceptance rate
cat("Acceptance rate:", mean(accept), "\n")

# Plot the samples and the target distribution
hist(x, breaks = 30, freq = FALSE, main = "mixture of two normal distributions")
curve(f(x, mu1, mu2, sigma1, sigma2, w), add = TRUE, col = 'red', lwd = 2)

library(mixtools)
fit <- normalmixEM(x, k = 2) #we are fitting a mixture of two Gaussian distributions
# Plot the fitted mixture model
curve(fit$lambda[1] * dnorm(x, fit$mu[1], sqrt(fit$sigma[1]^2)) +
        fit$lambda[2] * dnorm(x, fit$mu[2], sqrt(fit$sigma[2]^2)),
      from = min(x), to = max(x), add = TRUE, col = "blue")
