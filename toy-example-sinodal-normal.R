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
  # We check both endpoints simultaneously to avoid the need of two loops.
  while (target(x - a) >= target(x) * y || target(x + b) >= target(x) * y) {
    a <- 2 * a
    b <- 2 * b
  }
  # while (target(x + b) >= target(x) * y) {
  #   b <- 2 * b
  # }
  return(c(x - a, x + b))
}
# Doing a single slice sample.
SliceSampling <- function (x, alpha=1) {
  y <- runif(1)
  c <- RandomInterval(x, y, alpha)
  z <- runif(1, c[1], c[2])  # runif(1, -4, 4)
  while (target(z) < target(x) * y) {
    c <- RandomInterval(x, y, alpha)
    z <- runif(1, c[1], c[2])
  }
  return(z)
}

# Testing the MCMC method including a histogram and the real density.
# the MH algorithm and SliceSampling can simple be exchanged.
T <- 10^4
x <- rep(3.14, T)
for (t in 2:T) x[t] <- SliceSampling(x[t-1], 2)
# Throwing away the burn in period. Somehow this does make it worse...
# x <- x[1:T/10]
hist(x, breaks=seq(min(x), max(x), length=100), freq=FALSE, ylim=c(0, 0.7))
y <- seq(min(x), max(x), length=500)
z <- target(y) / Z
lines(y, z, col="red", lwd=2)

# Calculating the acceptence rate for the MH algorithm; 25% is desired
sum(x[-1] != x[1:(length(x) - 1)])/(T - 1)

# Multivariate setting.
# Define the unnormalised density.
target2 <- function(x){
  # sin(x[1])^2 * sin(2 * x[2])^2 * dnorm(x[1]) * dnorm(x[2])
  sin(x)^2 * sin(2 * x)^2 * dnorm(x)
}
# Load library for multidimensional integration to compute the normalisation
# constant.
library(cubature)
Z <- hcubature(target2, c(-100, -100), c(100, 100))[[1]]

# Implement the propose and reject step.
# Load library for multivariate normal.
library(MASS)
metropolis2 <- function(x, alpha=1){
  d <- length(x)
  y <- mvrnorm(1, x, diag(rep(alpha, d), d))
  if (runif(1) > target2(y) / target2(x)) y <- x
  return(y)
}

# Now we turn towards slice sampling.
# Proposing a random interval that includes the slice.
RandomInterval2 <- function (x, y, alpha=1) {
  d <- length(x)
  a <- rexp(1, rate=alpha)  # rexp(length(x), rate=alpha)
  b <- rexp(1, rate=alpha)  # rexp(length(x), rate=alpha)
  # We check both endpoints simultaneously to avoid the need of two loops.
  while (target2(x - rep(a, d)) >= target2(x) * y || target2(x + rep(b, d)) >= target2(x) * y) {
    a <- 2 * a
    b <- 2 * b
  }
  # while (target(x + b) >= target(x) * y) {
  #   b <- 2 * b
  # }
  return(matrix(c(x - a, x + b), d))
}

# Doing a single slice sample.
SliceSampling2 <- function (x, alpha=1) {
  d <- length(x)
  y <- runif(1)
  c <- RandomInterval2(x, y, alpha)
  z <- runif(d, c[, 1], c[, 2])  # runif(1, -4, 4)
  while (target2(z) < target2(x) * y) {
    c <- RandomInterval2(x, y, alpha)
    z <- runif(d, c[, 1], c[, 2])
  }
  return(z)
}

T <- 10^4
x <- matrix(rep(.2, T), 1)
for (t in 2:T) x[, t] <- metropolis2(x[, t-1], 2)
plot(t(x), pch=16, col='black', cex=0.5)
k <- kde2d(x[1, ], x[2, ], n=200, lims = c(-4, 4, -4, 4))
image(k, col=r, xlim=c(-4, 4), ylim=c(-4, 4))
y1 <- seq(min(x[2, ]), max(x[2, ]), length=500)
z1 <- sin(2 * y1)^2 * dnorm(y1) * 3
lines(z1 - rep(3.5, 500), y1, col="red", lwd=2)
y2 <- seq(min(x[1, ]), max(x[1, ]), length=500)
z2 <- sin(y2)^2 * dnorm(y2) * 3
lines(y2, z2 - rep(3.5, 500), col="red", lwd=2)


# hist(x, breaks=seq(min(x), max(x), length=100), freq=FALSE)


source("http://bioconductor.org/biocLite.R")
biocLite("hexbin")
y
install.packages("hexbin")
y
library(hexbin)
a
y
install.packages("gplots")
y
library(gplots)

library(MASS)

x <- rnorm(mean=1.5, 5000)
y <- rnorm(mean=1.6, 5000)
df <- data.frame(x,y)
plot(df, pch=16, col='black', cex=0.5)
k <- kde2d(df$x, df$y, n=200)
image(k, col=r)
hist2d(1)

library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

library(hexbin)
# Create hexbin object and plot
h <- hexbin(df)
plot(h)
plot(h, colramp=rf)
 