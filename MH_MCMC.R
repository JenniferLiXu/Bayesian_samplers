#Define target distribution(sinus function)
target <- function(x) {
  sin(x)
}

# Define the proposal distribution(a mixture of two gaussian distribution)
proposal <- function(x, mu, sigma, w) {
  if (runif(1) < w) {
    rnorm(1, mean = x, sd = 1)
  } else {
    rnorm(1, mean = x - mu, sd = sigma)
  }
}

likeli <- function(x, mu, m, sigma, w) {
  w*dnorm(x, mean = m, sd = 1) + (1-w)*dnorm(x-mu, mean = m, sd = sigma)
}

#initialize
set.seed(121) 
n <- 10000
x <- numeric(n)
accept <- numeric(n)

mu <- 2
sigma <- 0.5
w <- 0.7

#generate samples using the Metropolis-Hastings algorithm
for (i in 2:n) {
  # Generate candidate sample
  y <- proposal(x[i-1],  mu, sigma, w)
  
  # Calculate acceptance ratio
    alpha <- min(1, target(y) * likeli(y, mu, x[i-1],sigma, w) /
                 (target(x[i-1])*likeli(x[i-1], mu, y,sigma, w) ))

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

par(mfrow = c(1,1))
# Plot the samples and the target distribution
hist(x, breaks = 100, freq = FALSE, main = "mixture of two normal distributions")
curve(target(x), add = TRUE, col = 'red', lwd = 2)
abline(v = mean(x[-(1:burnIn)]))

library(mixtools)
fit <- normalmixEM(x, k = 2) #we are fitting a mixture of two Gaussian distributions
# Plot the fitted mixture model
curve(fit$lambda[1] * dnorm(x, fit$mu[1], sqrt(fit$sigma[1]^2)) +
        fit$lambda[2] * dnorm(x, fit$mu[2], sqrt(fit$sigma[2]^2)),
      from = min(x), to = max(x), add = TRUE, col = "blue")


#plot(x[-(1:burnIn)], type = "l", xlab="" , main = "Chain values of x", )



