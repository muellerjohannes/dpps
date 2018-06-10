# In this example we sample points on a (discrete) line according to a DPP
# using the dual sampling algorithm. We model L directly and via the
# quality-diversity decomposition. We plot and compare the patterns to
# uncorrelated points i.e. to a Poisson point process.

# Minimal example _____________________________________________________________
n <- 3
C <- matrix(c(2,1,0,1,2,0,0,0,2), nrow=n)

# Modelling phi and q _________________________________________________________
m <- 99
n <- m + 1
q <- rep(sqrt(m), n)
phi <- rep(0, n^2)
for (i in 1:n) {
  for (j in 1:n) {
    phi[(i - 1) * n + j] <- dnorm((i - j) / sqrt(m))
  }
}
phi <- matrix(phi, ncol=n)

# General part, define C ______________________________________________________
d <- length(phi) / n
for (i in 1:d) {
  phi[i, ] <- sum(phi[i, ]^2)^(-1/2) * phi[i, ]
}
B <- t(phi) * q
C <- t(B) %*% B

# Compute the eigendecomposition of C, set near zero eigenvalues to zero ______
edc <- eigen(C)
lambda <- edc$values
lambda[abs(lambda) < 10^(-9)] <- 0
mean <- sum(lambda / (1 + lambda))
eigenvectors <- edc$vectors

# Sample and plot things ______________________________________________________
time <- proc.time()
dataDPP <- sort(DualSamplingDPP(lambda, eigenvectors, C, B))
sort(DualSamplingDPP(lambda, eigenvectors, C, B))
pointsDPP <- t(CoordinatesNew(dataDPP, m))
plot(pointsDPP, xlim=0:1, ylim=0:1, xlab="", ylab="", xaxt='n', yaxt='n', asp=1)
proc.time() - time
pointsDPP <- DualSamplingDPP(lambda, eigenvectors, C, B)
plot(rep(1, length(pointsDPP)), pointsDPP,
     ylim=c(1, n), xlim=c(.4, 3.2), xaxt='n', ylab="Points", xlab="")

# Remove all objects apart from functions
rm(list = setdiff(ls(), lsf.str()))
