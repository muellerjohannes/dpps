# In this example we sample points on a (discrete) line according to a DPP
# We model L directly and via the quality-diversity decomposition. We plot and
# compare the patterns to uncorrelated points i.e. to a Poisson point process.

# Define the function that computes L from q and phi as well as the truncated
# eigendecomposition and a poisson points process of given cardinality.
DefineS <- function (q, phi) {
  n <- length(q)
  for (i in 1:n) {
    phi[, i] <- sum(phi[, i]^2)^(-1/2) * phi[, i]
  }
  S <- t(phi) %*% phi
  return(S)
}
DefineL <- function (q, phi) {
  n <- length(q)
  for (i in 1:n) {
    phi[, i] <- sum(phi[, i]^2)^(-1/2) * phi[, i]
  }
  S <- t(phi) %*% phi
  L <- t(q * S) * q
  return(L)
}
TruncatedEigendecomposition <- function (L) {
  edc <- eigen(L)
  lambda <- edc$values
  lambda[lambda < 10^(-9)] <- 0
  eigenvectors <- edc$vectors
  return(list(lambda, eigenvectors))
}
DefinePoissonPoints <- function (mean, n) {
  lambda <- rep(mean / n / (1 - mean / n), n)
  eigenvectors <- diag(rep(1, n))
  return(list(lambda, eigenvectors))
}

# Define the size of the ground set
m <- 99  # In case of the 0-1 sequence: 29
n <- m + 1

# Modelling phi and q _________________________________________________________
# Constant qualities
q <- rep(10^2, n)  # 0-1 sequences: rep(10^2, n)
# Log linear qualities
for (i in 1:n) {
  q[i] <- 10^2 * sqrt(m) * exp(-0.2 * abs(i - 50.5))
}
# Define the diversity features phi
phi <- rep(0, n^2)
for (i in 1:n) {
  for (j in 1:n) {
    phi[(i - 1) * n + j] <- dnorm((i - j), sd=10)  # 0-1 sequences: dnorm(0.437 * (i - j))
  }
}
phi <- matrix(phi, ncol=n)
# Define L
L <- DefineL(q, phi)
# Alternative direct definition of L
L <- rep(0, n^2)
for (i in 1:n) {
  for (j in 1:n) {
    L[(i - 1) * n + j] <- dnorm((i-j) * n^(-1/4))
  }
}
L <- matrix(L, nrow=n)

# Compute the eigendecomposition and set up the poisson points
edc <- TruncatedEigendecomposition(L)
lambda <- edc[1][[1]]
eigenvectors <- edc[2][[1]]
edc <- DefinePoissonPoints(sum(lambda / (1 + lambda)), n)
lambda2 <- edc[1][[1]]
eigenvectors2 <- edc[2][[1]]

# Sample and plot things ______________________________________________________
# Minimal example

# 0-1 sequences
x <- sort(SamplingDPP(lambda, eigenvectors))
as.integer(1:n %in% x)
y <- sort(SamplingDPP(lambda2, eigenvectors2))
as.integer(1:n %in% y)

# Sample from both point processes and plot the points on the line
pointsDPP <- SamplingDPP(lambda, eigenvectors)
pointsPoisson <- SamplingDPP(lambda2, eigenvectors2)
plot(rep(1, length(pointsDPP)), pointsDPP,
     ylim=c(1, n), xlim=c(.4, 3.2), xaxt='n', ylab="Points", xlab="")
points(rep(2, length(pointsPoisson)), pointsPoisson, pch=5)
legend("topright", inset=.05, legend=c("DPP", "Poisson"), pch=c(1, 5))

# Remove all objects apart from functions
rm(list = setdiff(ls(), lsf.str()))

# Set up dual sampling
# edc <- eigen(C)
# lambda <- edc$values
# lambda[abs(lambda) < 10^(-9)] <- 0
# mean <- sum(lambda / (1 + lambda))
# eigenvectors <- edc$vectors

# Dual sampling
# time <- proc.time()
# dataDPP <- sort(DualSamplingDPP(lambda, eigenvectors, C, B))
# sort(DualSamplingDPP(lambda, eigenvectors, C, B))
# pointsDPP <- t(CoordinatesNew(dataDPP, m))
# plot(pointsDPP, xlim=0:1, ylim=0:1, xlab="", ylab="", xaxt='n', yaxt='n', asp=1)
# proc.time() - time
# pointsDPP <- DualSamplingDPP(lambda, eigenvectors, C, B)
# plot(rep(1, length(pointsDPP)), pointsDPP,
#      ylim=c(1, n), xlim=c(.4, 3.2), xaxt='n', ylab="Points", xlab="")
