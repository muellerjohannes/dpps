# NEEDS: Sampling algorithm, declaration of the points in the square
# TODO: Maybe do the gradient descent directly over the representation
# od the gradient

# With this toy example we aim to perform the first learning of paramters
# associated to a kernel of a DPP. More precisely we will generate our own
# data of points on a line and estimate the standard deviation of the Gaussian
# kernels that introduce the diversifying behaviour.

# Generation of data
time <- proc.time()
T <- 5
data <- rep(list(0), T)
for (i in 1:T) {
  data[[i]] <- sort(SamplingDPP(lambda, eigenvectors))
}
proc.time() - time
# data

# Define the similarity S, L and the loss in dependency of the parameter sigma
Similarity <- function(sigma) {
  phi <- rep(0, n^2)
  for (i in 1:n) {
    for (j in 1:n) {
      phi[(i - 1) * n + j] <- dnorm((i - j), sd=sigma)  # 0-1 sequences: devide by 2
    }
  }
  phi <- matrix(phi, ncol=n)
  for (i in 1:n) {
    phi[, i] <- sum(phi[, i]^2)^(-1/2) * phi[, i]
  }
  S <- t(phi) %*% phi
  return(S)
}
LFunction <- function(sigma) {
  S <- Similarity(sigma)
  return(t(q * S) * q)
}
Loss <- function(sigma) {
  T <- length(data)
  # Sum this over all data entries!
  x <- 0
  for (i in 1:T) {
    A <- data[[i]]
    x <- x + log(det(matrix(Similarity(sigma)[A, A], length(A))))
  }
  return(- x + T * log(det(diag(rep(1, n)) + LFunction(sigma))))
}

# Parameter estimations
time <- proc.time()
sol <- nlm(Loss, 10)
proc.time() - time
sol$estimate
sol
warnings()
Loss(10)
Loss(7)

# Algorithm for the gradient of the loss function


# Remove all objects apart from functions
rm(list = setdiff(ls(), lsf.str()))
