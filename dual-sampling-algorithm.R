# TODO: Speed up selection of linear independent vectors
# TODO: How can I avoid double selections of elements?

ScalarProduct <- function(x, y, C) {
  return(t(x) %*% C %*% y)
}

DualSamplingDPP <- function (lambda, eigenvectors, B, C) {
  # First part of the algorithm, doing the selection of the eigenvectors
  D = length(lambda)
  J <- runif(D) <= lambda/(1 + lambda)
  k <- sum(J)
  # V <- matrix(rep(0, D * k))
  # for (i in 1:D) {
  #   V[, i] <- eigenvectors[]
  # }
  V <- matrix(t(t(eigenvectors[, J]) / diag(ScalarProduct(eigenvectors[, J], eigenvectors[, J], C))), nrow=D)
  Y <- rep(0, k)
  
  # Second part of the algorithm, the big while loop
  while (k > 0) {
    # Calculating the weights and selecting an item i according to them
    wghts <- k^(-1) * colSums((t(V) %*% B)^2)
    i <- sample(D, 1, prob=wghts)
    Y[k] <- i
    if (k == 1) break
    
    # Projecting B_i onto the span of V
    help <- V %*% t(V) %*% B[, i]
    help <- sum(help^2)^(-1/2) * help
    
    # Projecting the elements of V onto the subspace orthogonal to B_i/the projection help of B_i
    V <- V - help %*% t(t(V) %*% help)
    
    # Orthonormalize V and set near zero entries to zero
    # V[abs(V) < 10^(-9)] <- 0  # This might delete all entries, if C is big!
    j <- 1
    while(j <= k) {
      help2 <- rep(0, D)
      m <- 1
      while (m <= j - 1) {
        help2 <- help2 + ScalarProduct(V[, j], V[, m], C) * V[, m]
        m <- m + 1
      }
      # help2 <- (matrix(V[, 1:(j - 1)], nrow=N) %*% colSums(matrix(V[, j] * V[, 1:(j - 1)], nrow=N)))[, 1]
      V[, j] <- V[, j] - help2
      if (ScalarProduct(V[, j], V[, j], C) > 0) {
        V[, j] <- ScalarProduct(V[, j], V[, j], C)^(1 / 2) * V[, j]
      }
      j <- j + 1
    }
    # V[abs(V) < 10^(-9)] <- 0
    
    # Selecting a linear independent set in V
    k <- k - 1
    q <- qr(V)
    V <- matrix(V[, q$pivot[seq(k)]], ncol=k)
  }
  return(Y)
}