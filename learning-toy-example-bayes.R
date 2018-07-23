# NEEDS: MCMC algorithm, learning toy example

# We want to introduce the first example of Bayesian parameter estimation for
# DPPs. We start by estimating the log linearity constant of the qualities - a
# parameter for which we've already successfully done MLE.

# Define the unnormalised density.
target <- function (theta) {
  T <- length(data)
  x <- 1
  a <- det(diag(rep(1, n)) + LFunction(theta)) / 10^27
  for (t in 1:T) {
    A <- data[[t]]
    x <- x * exp(2 * sum(theta * Feature(A))) * det(S[A, A]) / a
  }
  x <- x
  return(x)
}

# Run MCMC and create a heatmap that also shows the MLE.
time <- proc.time()
x <- MCMC(target, c(-8, 2), MH=FALSE, 10^3, alpha=3)  # 3 or 0.1
proc.time() - time
plot(t(x), pch=16, col='black', cex=0.5, xlim=c(-9.5, -6.5), ylim=c(1.5, 2.8))
k <- kde2d(x[1, ], x[2, ], n=1000, lims = c(-11, -3, -1, 3))
image(k, col=r, xlim=c(-9.5, -6.3), ylim=c(1.4, 2.8))
points(z[1], z[2], pch=4, lwd=3, col="green")
lines(c(z[1], z[1]), c(1, 3), type="l", col="red", lwd="1.5")
lines(c(-10, -6), c(z[2], z[2]), type="l", col="red", lwd="1.5")

# Calculating the acceptence rate for the MH algorithm; around 25% is desired
sum(x[, -1] != x[, 1:(length(x)/2 - 1)])/(length(x) - 1)
