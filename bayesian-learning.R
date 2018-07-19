# The goal is to implement the Bayesian approach to the kernel estimation og
# discrete DPPs presented.

# Define the elementary kernel L
L <- matrix(c(2, 1, 1, 2), 2)
L

# Function for the unnormalised posterior
UnnormalisedPosterior <- function (theta, data) {
  T <- length(data)
  x <- 1
  for (i in 1:T) {
    x <- x * det(theta[data[[i]], data[[i]]])
  }
  return(x)
}

L
T <- 3
data <- rep(list(0), T)
for (i in 1:T) {
  data[[i]] <- sort(SamplingDPP(lambda, eigenvectors))
}
data
length(data)
data[[1]]
L[data[[1]], data[[1]]]

