# TODO: Implement the methods in the multivariate case.
# Add the burn in period for the MH random walk. Also adjust the proposal.
# Creating a first toy example for MCMC methods.
# First we implement the Metropolis-Hastings algorithm.

# Computation of the actual partition function.
target <- function(x){
  sin(x)^2 * sin(2 * x)^2 * dnorm(x)
}
Z <- integrate(target, -100, 100)[[1]]

# Implement the propose and reject step.
metropolis <- function(x, alpha=1){
  y <- rnorm(1, mean=x, sd=alpha)  # runif(1, x - alpha, x + alpha)
  if (runif(1) > target(y) / target(x)) y <- x
  return(y)
}

# Now we turn towards slice sampling.
# Proposing a random interval that includes the slice.
RandomInterval <- function (x, y, alpha=1) {
  a <- rexp(1, rate=alpha)
  b <- rexp(1, rate=alpha)
  while (target(x - a) >= target(x) * y) {
    a <- 2 * a
  }
  while (target(x + b) >= target(x) * y) {
    b <- 2 * b
  }
  return(c(x - a, x + b))
}
# Doing a single slice sample.
SliceSampling <- function (x, alpha=1) {
  y <- runif(1)
  c <- RandomInterval(x, y, alpha)
  z <- runif(1, c[1], c[2])  # runif(1, -4, 4)
  while (target(z) < target(x) * y) {
    z <- runif(1, -4, 4)
  }
  return(z)
}

# Testing the MCMC method including a histogram and the real density.
# the MH algorithm and SliceSampling can simple be exchanged.
T <- 10^4
x <- rep(3.14, T)
for (t in 2:T) x[t] <- SliceSampling(x[t-1], 2)
hist(x, breaks=seq(min(x), max(x), length=100), freq=FALSE)
y <- seq(min(x), max(x), length=500)
z <- target(y) / Z
lines(y, z, col="red", lwd=2)

# Calculating the acceptence rate for the MH algorithm; 25% is desired
sum(x[-1] != x[1:(length(x) - 1)])/(T - 1)