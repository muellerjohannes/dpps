# TODO: Add the burn in period for the MH random walk. Understand why this
# does not improve the accuracy. Also adjust the proposal.

# README: Creating some toy examples for MCMC methods.

# First we implement the Metropolis-Hastings algorithm. We implement the
# propose and reject step. We use a Gaussian as a proposal with covariance
# matrix alpha times the identity.
# Load library for multivariate normal.
library(MASS)
Metropolis <- function(x, f, alpha=1){
  d <- length(x)
  if (length(alpha) == 1) {
    alpha <- diag(rep(alpha, d), d)
  }
  y <- mvrnorm(1, x, alpha)
  if (runif(1) * f(x) > f(y)) y <- x
  return(y)
}

# Now we turn towards slice sampling. Proposing a random interval that includes
# the slice. We use an exponential random variable to define the width of the
# interval.
RandomInterval <- function (x, y, f, alpha=1) {
  c <- f(x)
  # We make the interval the same length in every dimension.
  a <- rexp(1, rate=alpha)  # rexp(length(x), rate=alpha)
  b <- rexp(1, rate=alpha)  # rexp(length(x), rate=alpha)
  # One can check both endpoints simultaneously to avoid the need of two loops.
  while (TRUE) {
    help <- f(x - a)
    if (is.nan(help) || help < c * y) break
    a <- 2 * a
  }
  while (TRUE) {
    help <- f(x + b)
    if (is.nan(help) || help < c * y) break
    b <- 2 * b
  }
  return(matrix(c(x - a, x + b), length(x)))
}
# Doing a single slice sample.
SliceSampling <- function (x, f, alpha=1) {
  d <- length(x)
  a <- f(x)
  y <- runif(1)
  c <- RandomInterval(x, y, f, alpha)
  z <- runif(d, c[, 1], c[, 2])  # runif(1, -4, 4)
  while (TRUE) {
    help <- f(z)
    if (is.nan(help) || help < a * y) z <- runif(d, c[, 1], c[, 2])
    else break
    # c <- RandomInterval(x, y, f, alpha)
    # z <- runif(d, c[, 1], c[, 2])
    # z <- runif(d, c(-10, 0), c(-5, 4))
  }
  return(z)
}

# Implementing the MCMC method. The function needs the unnormalised density f,
# a starting value x0, sample size T whether it should be MH or Slice Sampling
# and the parameter alpha, which either specifies the variance of the proposal
# which is multivariate normal or the rate of the exponential random variable
# which defines the thickness of the random interval.
MCMC <- function (f, x0, T=10^3, MH=TRUE, alpha=1) {
  d <- length(x0)
  x <- matrix(rep(x0, T), d)
  if (MH) {
    for (t in 2:T) x[, t] <- Metropolis(x[, t-1], f, alpha) 
  }
  else {
    # Check whether starting value is impossible. In this case the slice is the
    # whole space and hence the endpoints of the random interval will diverge.
    while (is.nan(f(x0)) || f(x0)==0) {
      x0 <- mvrnorm(1, x0, diag(rep(alpha, d), d))
    }
    x[, 1] <- x0
    for (t in 2:T) x[, t] <- SliceSampling(x[, t-1], f, alpha)
  }
  return(x)
}

# One dimensional case, making nice pictures also.
# Define the unnormalised density we want to sample from.
target <- function(x){
  sin(x)^2 * sin(2 * x)^2 * dnorm(x)
}
# Load library for multidimensional integration to compute the normalisation
# constant.
library(cubature)
Z <- hcubature(target, -20, 20)[[1]]

x <- MCMC(target, 1, MH=TRUE, 10^4, alpha=1)
hist(x, breaks=seq(min(x), max(x), length=100), freq=FALSE, ylim=c(0, .7))
y <- seq(min(x), max(x), length=500)
z <- target(y) / Z
lines(y, z, col="red", lwd=2)

# Calculating the acceptence rate for the MH algorithm; around 25% is desired
sum(x[-1] != x[1:(length(x) - 1)])/(length(x) - 1)

# Two dimensional toy example with a similar density. The points are plotted
# and a heat map is created which also shows the marginal densities.
target2 <- function(x){
  sin(x[1])^2 * sin(2 * x[2])^2 * dnorm(x[1]) * dnorm(x[2])
}
x <- MCMC(target2, c(1, 1), MH=TRUE, 10^4, alpha=3)
plot(t(x), pch=16, col='black', cex=0.5)
# Color housekeeping and loading library for heat map.
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)
library(MASS)
k <- kde2d(x[1, ], x[2, ], n=200, lims = c(-4.3, 4, -4, 4))
image(k, col=r, xlim=c(-4.3, 4), ylim=c(-4, 4))
# Adding the marginal densities.
y1 <- seq(-4, 4, length=500)
z1 <- sin(2 * y1)^2 * dnorm(y1) * 3
lines(z1 - rep(3.8, 500), y1, col="red", lwd=2)
y2 <- seq(-4.3, 4, length=500)
z2 <- sin(y2)^2 * dnorm(y2) * 3
lines(y2, z2 - rep(3.5, 500), col="red", lwd=2)

# Remove all objects apart from functions
rm(list = setdiff(ls(), lsf.str()))
rm(list = ls())
