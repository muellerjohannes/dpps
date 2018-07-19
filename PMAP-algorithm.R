# The goal is to implement the PMAP algorithm, that reconstructs the
# determinantal equivalence class of a symmetric matrix given its principle
# minors.

K <- matrix(c(2, -1, 0, -1, 2, 1, 0, 1, 2), 3)/4
N <- sqrt(length(K))
K

# Construction the estimate up to the signs of the off diagonal values
A <- diag(diag(K))
for (i in 1:(N - 1)) {
  for (j in (i + 1):N) {
    A[[i, j]] <- sqrt(A[[i, i]] * A[[j, j]] - det(K[c(i, j), c(i, j)]))
    A[[j, i]] <- A[[i, j]]
  }
}

A
