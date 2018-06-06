# Different kernels L _________________________________________________________
# _____________________________________________________________________________

# Minimal example _____________________________________________________________
n <- 3
L <- matrix(c(2,1,0,1,2,0,0,0,2), nrow=n)

# Points on a line ____________________________________________________________
n <- 100
L <- rep(0, n^2)
for (i in 1:n) {
  for (j in 1:n) {
    L[(i - 1) * n + j] <- dnorm((i-j) * n^(-1/4))
  }
}
L <- matrix(L, nrow=n)

# Points in the square ________________________________________________________
# Define the coordinates of a point
Coordinates <- function(i, n) {
  y1 <- floor((i - 1) / (n + 1))
  x1 <- i - 1 - (n + 1) * y1
  return (c(x1, y1)/n)
}
Distance <- function (i, j, n) {
  return (sqrt(sum((Coordinates(i, n) - Coordinates(j, n))^2)))
}
m <- 19
n <- (m + 1)^2
L <- rep(0, n^2)
for (i in 1:n) {
  for (j in 1:n) {
    L[(i - 1) * n + j] = n^2 * dnorm(Distance(i, j, m))
  }
}
L <- matrix(L, nrow=n)

# Modelling phi and q _________________________________________________________
# Points on the line.
m <- 99
n <- m + 1
q <- rep(0, n)
for (i in 1:n) {
  q[i] <- sqrt(m)   # 10^2 * sqrt(m) * exp(-0.3 * abs(i - 50.5))
}
phi <- rep(0, n^2)
for (i in 1:n) {
  for (j in 1:n) {
    phi[(i - 1) * n + j] <- dnorm((i - j) / sqrt(m))
  }
}

# Points in the square.
m <- 29
n <- (m + 1)^2
q <- rep(0, n)
for (i in 1:n) {
  q[i] <- sqrt(m)
}
phi <- rep(0, n^2)
time <- proc.time()
for (i in 1:n) {
  for (j in 1:n) {
    phi[(i - 1) * n + j] <- dnorm(sqrt(m) * Distance(i, j, m))
  }
}
proc.time() - time

# Quality diversity decomposition with small D
phi <- rep(0, n * sqrt(n) * 10)
time <- proc.time()
for (i in 1:n) {
  for (j in 1:(sqrt(n) * 10)) {
    phi[(i - 1) * sqrt(n) + j] <- dnorm(sqrt(m) *
                                          Distance(i, sqrt(n) / 10 * j, m))
  }
}
proc.time() - time
phi <- matrix(phi, ncol=n)

# Log linear quality for the points on the line
m <- 99
n <- m + 1
q <- rep(0, n)
for (i in 1:n) {
  q[i] <- 10^2 * sqrt(m) * exp(-0.2 * abs(i - 50.5))
}
phi <- rep(0, n^2)
for (i in 1:n) {
  for (j in 1:n) {
    phi[(i - 1) * n + j] <- dnorm(2 * (i - j) / sqrt(m))
  }
}

# General part.
phi <- matrix(phi, nrow=n)
for (i in 1:n) {
  phi[, i] <- sum(phi[, i]^2)^(-1/2) * phi[, i]
}
L <- rep(0, n^2)
time <- proc.time()
for (i in 1:n) {
  for (j in 1:n) {
    L[(i - 1) * n + j] <- q[i] * q[j] * sum(phi[, i] * phi[, j])
  }
}
proc.time() - time
L <- matrix(L, nrow=n)

# Compute the eigendecomposition, set near zero eigenvalues to zero and
# set up poisson point process with same expected cardinality _________________
edc <- eigen(L)
lambda <- edc$values
lambda[abs(lambda) < 10^(-9)] <- 0
min(lambda)
mean <- sum(lambda / (1 + lambda))
eigenvectors <- edc$vectors
D <- diag(rep(mean / n / (1 - mean / n), n))
edc2 <- eigen(D)
lambda2 <- edc2$values
sum(lambda2 / (1 + lambda2))
eigenvectors2 <- edc2$vectors

# Sample and plot things ______________________________________________________
# _____________________________________________________________________________

# Minimal example
SamplingDPP(lambda, eigenvectors)

# Sample from both point processes and plot the points on the line
pointsDPP <- SamplingDPP(lambda, eigenvectors)
pointsPoisson <- SamplingDPP(lambda2, eigenvectors2)
plot(rep(1, length(pointsDPP)), pointsDPP,
     ylim=c(1, n), xlim=c(.4, 3.2), xaxt='n', ylab="Points", xlab="")
points(rep(2, length(pointsPoisson)), pointsPoisson, pch=5)
legend("topright", inset=.05, legend=c("DPP", "Poisson"), pch=c(1, 5))

# Sample from both point processes and plot the points in the square
# par(mfrow = c(1,1))
dataDPP <- sort(SamplingDPP(lambda, eigenvectors))
pointsDPP <- matrix(Coordinates(dataDPP, m), ncol=2)
plot(pointsDPP, xlim=0:1, ylim=0:1, xlab="", ylab="", asp=1)
dataPoisson <- sort(SamplingDPP(lambda2, eigenvectors2))
pointsPoisson <- matrix(Coordinates(dataPoisson, m), ncol=2)
plot(pointsPoisson, xlim=0:1, ylim=0:1,
     # xaxt='n', yaxt='n',
     xlab="", ylab="", asp=1)

# Remove all objects apart from functions
rm(list = setdiff(ls(), lsf.str()))
