log_euclidean_mean <- function(sigmaList) {
  n <- length(sigmaList)
  sigma <- diag(0, 3L)
  for (i in 1:n) {
    eig <- eigen(sigmaList[[i]], symmetric = TRUE)
    val <- log(eig$values)
    tmp <- eig$vectors %*% diag(val) %*% t(eig$vectors)
    sigma <- sigma + tmp
  }
  eig <- eigen(sigma / n, symmetric = TRUE)
  sigma <- eig$vectors %*% diag(exp(eig$values)) %*% t(eig$vectors)
  sigma
}

log_euclidean_distance <- function(sigma1, sigma2) {
  eig1 <- eigen(sigma1, symmetric = TRUE)
  eig2 <- eigen(sigma2, symmetric = TRUE)
  val1 <- log(eig1$values)
  val2 <- log(eig2$values)
  log1 <- eig1$vectors %*% diag(val1) %*% t(eig1$vectors)
  log2 <- eig2$vectors %*% diag(val2) %*% t(eig2$vectors)
  diff <- log1 - log2
  sum(diff^2)
}
