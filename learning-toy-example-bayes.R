# NEEDS: MCMC algorithm, learning toy example

# We want to introduce the first example of Bayesian parameter estimation for
# DPPs. We start by estimating the log linearity constant of the qualities - a
# parameter for which we've already successfully done MLE.

# Define the unnormalised density.
target <- function (theta) {
  T <- length(data)
  x <- 1
  a <- det(diag(rep(1, n)) + LFunction(theta))
  for (t in 1:T) {
    A <- data[[t]]
    x <- x * exp(2 * sum(theta * Feature(A))) * det(S[A, A]) / a
  }
  x <- x
  return(x)
}

# Run MCMC and create a heatmap that also shows the MLE.
time <- proc.time()
x <- MCMC(target, mle, MH=TRUE, 10^3, alpha=0.1)  # 3 or 0.1
proc.time() - time
plot(t(x2), pch=16, col='black', cex=0.5, xlim=c(-10, -7), ylim=c(1.4, 2.8))
lines(c(mle[1], mle[1]), c(-1, 3), type="l", col="red", lwd="1.5")
lines(c(-11, 1), c(mle[2], mle[2]), type="l", col="red", lwd="1.5")

# Calculating the acceptence rate for the MH algorithm; around 25% is desired
sum(x[, -1] != x[, 1:(length(x)/2 - 1)])/(length(x) - 2)

# Doing PCA with the first 100 samples in order to tune the proposal.
library(stats)
pc <- prcomp(t(x))
sd <- pc[[2]] %*% diag(c(pc[[1]][[1]], pc[[1]][[2]])^2) %*% pc[[2]]
time <- proc.time()
xnew <- MCMC(target, c(-8, 2), MH=TRUE, 10^4, alpha=sd)
proc.time() - time
k <- kde2d(xnew[1, ], xnew[2, ], n=500, lims = c(-11, -3, -1, 3))
image(k, col=r, xlim=c(-10, -6.8), ylim=c(1.4, 2.8))
points(mle[1], mle[2], pch=4, lwd=3, col="green")
