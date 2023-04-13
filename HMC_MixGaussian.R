hmc <- function(U, epsilon, L, current_q) {
  #input: 
  # U: the potential energy function 
  # epsilon: the step size 
  # L: the number of steps
  # current_q: the current position 
  chain <- numeric(n_iter)
  chain[1] <- current_q
  accept <- numeric(n_iter)
  
  for (i in 2:n_iter){
    q <- chain[i-1]
    p <- matrix(rnorm(length(q),0,1), ncol = 1) 
    current_p <- p
    
    # Leapfrog method: 
    # First make a half step for momentum 
    p <- p - epsilon / 2 * U(q, returnGrad = TRUE)  
    for (j in 1:L) {
      # Make a full step for the position
      q <- q + epsilon * p 
      # Make a full step for the momentum, except at end of trajectory
      # Make a half step for momentum at the end
      p <- p - epsilon * U(q, returnGrad = TRUE)/(1 + (j==L))
    }
  
    # Negate momentum at end of trajectory to make the proposal symmetric
    #p <- -p
    # Evaluate potential and kinetic energies at start and end of trajectory
    current_U <- U(chain[i-1], returnGrad = FALSE)
    current_K <- sum(current_p^2) / 2
    proposed_U <- U(q, returnGrad = FALSE)
    proposed_K <- sum(p^2) / 2
    
    if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K)) {
      chain[i] <- q
      accept[i] <- 1  #Accept
    } else {
      chain[i] <- chain[i-1] 
      accept[i] <- 0 #Reject
    }
  }
  return(list(chain = chain, accept_rate = accept))
  
}