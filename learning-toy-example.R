# With this toy example we aim to perform the first learning of paramters
# associated to a kernel of a DPP. More precisely we will generate our own
# data of points on a two dimensional grid with a log linear quality model
# and aim to estimate the log linearity parameter.

# Generation of data
T <- 10
data <- rep(list(0), T)
for (i in 1:T) {
  data[[i]] <- sort(SamplingDPP(lambda, eigenvectors))
}
# data

# Define the quality q, L, the feature sum and the loss in dependency of the
# parameter theta
# S <- t(phi) %*% phi
Quality <- function(theta) {
  return(exp(theta * DistanceNew(rep(5, n), 1:n, 2, m)))
}
LFunction <- function(theta) {
  return(t(t(Quality(theta) * S) * Quality(theta)))
} 
Feature <- function(A) {
  return(sum(DistanceNew(rep(5, length(A)), A, 2, m)))
}
Loss <- function(theta) {
  T <- length(data)
  # Sum this over all data entries!
  x <- 0
  for (i in 1:T) {
    A <- data[[i]]
    x <- x + log(det(matrix(LFunction(theta)[A, A], length(A))))
      # sum(theta * Feature(A)) + log(det(matrix(S[A, A], length(A))))
  }
  return(- x + T * log(det(diag(rep(1, n)) + LFunction(theta))))
}

# Parameter estimations
nlm(Loss, -4)

# Algorithm for the gradient of the loss function


# Remove all objects apart from functions
rm(list = setdiff(ls(), lsf.str()))
