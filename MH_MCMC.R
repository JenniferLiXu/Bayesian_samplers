#Define a mixture of two Gaussian distributions
f <- function(x, mu, sigma1, sigma2, w) {
  w * dnorm(x, mean = 2, sd = sigma1) + (1 - w) * dnorm(x, mean = 2 - mu, sd = sigma2)
}

#proposal distribution function
proposal <- function(x, sd) {
  rnorm(1, mean = x, sd = sd)
}

#initialize
set.seed(123)
n <- 10000
mu <- 2
sigma1 <- 1
sigma2 <- 2
w <- 0.6
sd <- 0.5

x <- numeric(n)
accept <- numeric(n)
f_value <- numeric(n)
f_value[1] <- f(0, mu, sigma1, sigma2, w)

#generate samples using the Metropolis-Hastings algorithm
for (i in 2:n) {
  # Generate candidate sample
  y <- proposal(x[i - 1], sd)
  f_proposed <- f(y, mu, sigma1, sigma2, w)
  # Calculate acceptance ratio
  alpha <- min(1, f_proposed / f_value[i - 1])

  # Decide whether to accept or reject the candidate sample
  if (runif(1) < alpha) {
    x[i] <- y
    accept[i] <- 1
    f_value[i] <- f_proposed
  }else {
    x[i] <- x[i - 1]
    accept[i] <- 0
    f_value[i] <- f_value[i - 1]
  }
}

#target \int \sin(x) p(x) dx
esti_sin <- mean(sin(x))

#numerical integral
num_sin <- integrate(function(x) sin(x) * f(x, mu, sigma1, sigma2, w), -Inf, Inf)

#difference between estimated value and numerical integral
abs(esti_sin - num_sin$value)

# Print the acceptance rate
cat("Acceptance rate:", mean(accept), "\n")

# Plot the samples and the target distribution
hist(x, breaks = 30, freq = FALSE, main = "mixture of two normal distributions")
curve(f(x, mu, sigma1, sigma2, w), add = TRUE, col = "red", lwd = 2)

library(mixtools)
#we are fitting a mixture of two Gaussian distributions
fit <- normalmixEM(x, k = 2)
# Plot the fitted mixture model
curve(fit$lambda[1] * dnorm(x, fit$mu[1], sqrt(fit$sigma[1]^2)) +
        fit$lambda[2] * dnorm(x, fit$mu[2], sqrt(fit$sigma[2]^2)),
      from = min(x), to = max(x), add = TRUE, col = "blue")

burnin <- 2000
plot(x[-(1:burnin)], type = "l", xlab = "", main = "Chain values of x", )
