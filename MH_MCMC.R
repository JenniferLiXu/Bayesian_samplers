MH_MCMC <- function(target, proposal, x_init, n_iter) {
  #target: the target distribution
  #proposal: the proposal distribution function
  #x_init: starting value of the chain
  #n_iter: number of iterations
  chain <- numeric(n_iter)
  accept <- numeric(n_iter)
  chain[1] <- x_init
  f_value <- target(x_init)$distri
  f_max <- f_value
  
  # Initialize the swap count and lengths of the chains between swaps
  swap_count <- 0
  chain_lengths <- integer(n_iter)
  chain_length <- 0
  
  for(i in 2:n_iter){
    # Generate candidate sample
    y <- proposal(chain[i - 1])
    f_proposed <- target(y)$distri
    # Calculate acceptance ratio
    alpha <- min(1, f_proposed / f_value)
    
    # Decide whether to accept or reject the candidate sample
    if (runif(1) < alpha) {
      chain[i] <- y
      accept[i] <- 1
      swap_occurred <- TRUE
      swap_count <- swap_count + 1
      f_value <- f_proposed
      f_max <- max(f_max, f_proposed)
    }else {
      swap_occurred <- FALSE
      chain[i] <- chain[i - 1]
      accept[i] <- 0
    }
    
    # Update the chain lengths
    if (swap_occurred) {
      chain_lengths[swap_count] <- chain_length
      chain_length <- 0
    } else {
      chain_length <- chain_length + 1
    }
  }
  return(list(chain = chain, f_max = f_max))
}

