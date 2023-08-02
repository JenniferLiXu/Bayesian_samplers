hmc <- function(U, epsilon, L, current_q) {
  #input: 
  # U: the potential energy function 
  # epsilon: the step size 
  # L: the number of steps
  # current_q: the current position 
  
  chain <- numeric()
  chain[1] <- current_q
  min_q <- current_q
  accept <- numeric(n_iter)
  swaps <- 0
  chain_lengths <- integer(n_iter - 1)
  min_U <- Inf  # Initialize maximum U
  intermediate_results <- numeric()
  
  for (i in 2:n_iter){
    q <- chain[i-1]
    p <- rnorm(length(q),0,1)
    current_p <- p
    
    # Leapfrog method: 
    # First make a half step for momentum 
    U_q <- U(q)
    p <- p - epsilon / 2 * U_q$grad
    for (j in 1:L) {
      # Make a full step for the position
      q <- q + epsilon * p 
      # Make a full step for the momentum
      # Make a half step for momentum at the end
      U_q_new <- U(q)
      p <- p - epsilon * U_q_new$grad/(1 + (j==L))
    }

    # Negate momentum at end of trajectory to make the proposal symmetric
    p <- -p
    # Evaluate potential and kinetic energies at start and end of trajectory
    #U_chain <- U(chain[i-1])
    current_U <- U_q$distri
    current_K <- sum(current_p^2) / 2
    proposed_U <- U_q_new$distri
    proposed_K <- sum(p^2) / 2
    
    accept_prob <- exp((current_U + current_K) - (proposed_U + proposed_K))
    
    if (!is.na(accept_prob) && !is.nan(accept_prob) && runif(1) < accept_prob) {
      #print(q)
      chain[i] <- q
      accept[i] <- 1  #Accept
      if (proposed_U < min_U) {
        min_q <- q
        min_U <- proposed_U # update max_q using max_U
      }
    } else {
      chain[i] <- chain[i-1] 
      accept[i] <- 0   #Reject
    }
    f_max <- target(min_q)$distri
  }
  
  return(list(chain = chain, accept_rate = accept, min_q = min_q, f_max = f_max))
  
}
