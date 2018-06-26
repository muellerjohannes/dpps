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
Quality <- function(theta) {
  return(rep(exp(theta), n))
  # return(exp(theta[1] * DistanceNew(rep(5, n), 1:n, 2, m) + theta[2]))
}
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
LFunction <- function(x) {
  theta <- x[1]
  sigma <- x[2]
  S <- Similarity(sigma)
  q <- Quality(theta)
  return(t(q * S) * q)
}
Loss <- function(x) {
  T <- length(data)
  # Sum this over all data entries!
  y <- 0
  for (i in 1:T) {
    A <- data[[i]]
    y <- y + log(det(matrix(LFunction(x)[A, A], length(A))))
  }
  return(- y + T * log(det(diag(rep(1, n)) + LFunction(x))))
}

# Parameter estimations
time <- proc.time()
sol <- nlm(Loss, c(2.3, 10))
proc.time() - time
sol$estimate
sol
warnings()
Loss(c(3, 10))
LFunction(2, 10)
A <- data[[1]]
log(det(matrix(LFunction(c(2, 10))[A, A], length(A))))

# Algorithm for the gradient of the loss function


# Remove all objects apart from functions
rm(list = setdiff(ls(), lsf.str()))
