MH_MCMC <- function(target, proposal, x_init, n_iter) {
  #target: the target distribution
  #proposal: the proposal distribution function
  #x_init: starting value of the chain
  #n_iter: number of iterations
  chain <- numeric(n_iter)
  accept <- numeric(n_iter)
  chain[1] <- x_init
  f_value <- target(x_init)
  
  for(i in 2:n_iter){
    # Generate candidate sample
    y <- proposal(chain[i - 1])
    f_proposed <- target(y, returnGrad = FALSE)
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