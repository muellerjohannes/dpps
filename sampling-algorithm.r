SamplingDPP <- function (lambda, eigenvectors) {
  # First part of the algorithm, doing the selection of the eigenvectors
  N = length(lambda)
  J <- runif(N) <= lambda/(1 + lambda)
  k <- sum(J)
  V <- matrix(eigenvectors[, J], nrow=N)
  Y <- c()
  
  # Second part of the algorithm, the big while loop
  while (k > 0) {
    # Calculating the weights and selecting an item i according to them
    wghts <- c()
    for (j in 1:N) {
      wghts <- append(wghts, k^(-1) * sum(V[j, ]^2))
    }
    # rm(j)
    i <- sample(N, 1, prob=wghts)
    # rm(wghts)
    Y <- append(Y, i)
    if (k == 1) break
    
    # Projecting e_i onto the span of V
    help <- rep(0, N)
    for (j in 1:k) {
      help <- help + V[i, j]*V[, j]
    }
    # rm(j)
    help <- sum(help^2)^(-1/2) * help
    
    # Projecting the elements of V onto the subspace orthogonal to e_i/the projection help of e_i
    for (j in 1:k) {
      V[, j] <- V[, j]-sum(V[, j] * help) * help
    }
    # rm(j)
    
    # Orthonormalize V and set near zero entries to zero
    V[abs(V) < 10^(-9)] <- 0
    j <- 1
    while(j <= k) {
      help2 <- rep(0, N)
      m <- 1
      while (m <= j - 1) {
        help2 <- help2 + sum(V[, j] * V[, m]) * V[, m]
        m <- m + 1
      }
      V[, j] <- V[, j] - help2
      if (sum(V[, j]^2) > 0) {
        V[, j] <- sum(V[, j]^2)^(-1/2) * V[, j]
      }
      # rm(m)
      j <- j + 1
    }
    V[abs(V) < 10^(-9)] <- 0
    
    # Selecting a linear independent set in V
    k <- k - 1
    q <- qr(V)
    V <- matrix(V[, q$pivot[seq(k)]], ncol=k)
  }
  # Remove variables
  # rm(list = setdiff(ls(), lsf.str()))  # Is this necessary? I think they are removed automatically
  return(Y)
}
