#' Addition in Bayes space
#'
#' @param isigma1 A \code{3x3} symmetric positive definite matrix.
#' @param isigma2 A \code{3x3} symmetric positive definite matrix.
#' @param neutral_variance The neutral element of the Bayes space (default: 1e5).
#'
#' @return A \code{3x3} symmetric positive definite matrix.
#' @export
#'
#' @examples
#' D1 <- diag(c(1.7, 0.3, 0.1)) * 1e-3
#' D2 <- diag(c(0.3, 1.7, 0.1)) * 1e-3
#' bayes_add(D1, D2)
bayes_add <- function(isigma1, isigma2, neutral_variance = neutral_element) {
  isigma1 + isigma2 - diag(1 / neutral_variance, 3L)
}

#' Scalar multiplication in Bayes space
#'
#' @param isigma A \code{3x3} symmetric positive definite matrix.
#' @param alpha A scalar.
#' @param neutral_variance The neutral element of the Bayes space (default: 1e5).
#'
#' @return A \code{3x3} symmetric positive definite matrix.
#' @export
#'
#' @examples
#' D <- diag(c(1.7, 0.3, 0.1)) * 1e-3
#' bayes_multiply(D, 2)
#' bayes_add(D, D)
bayes_multiply <- function(isigma, alpha, neutral_variance = neutral_element) {
  alpha * isigma - (alpha - 1) * diag(1 / neutral_variance, 3L)
}

#' Mean element in Bayes space
#'
#' @param isigmaList A list of \code{3x3} symmetric positive definite matrices.
#' @param neutral_variance The neutral element of the Bayes space (default: 1e5).
#'
#' @return A \code{3x3} symmetric positive definite matrix.
#' @export
#'
#' @examples
#' D <- list(D1 = diag(c(1.7, 0.3, 0.1)) * 1e-3,
#'           D2 = diag(c(0.3, 1.7, 0.1)) * 1e-3)
#' bayes_mean(D)
bayes_mean <- function(isigmaList, neutral_variance = neutral_element) {
  n <- length(isigmaList)
  w <- 1 / n
  isigma <- diag(1 / neutral_variance, 3L)
  for (i in 1:n) {
    tmp <- bayes_multiply(isigmaList[[i]], w, neutral_variance)
    isigma <- bayes_add(isigma, tmp, neutral_variance)
  }
  isigma
}

#' Distance in Bayes space
#'
#' @param isigma1 A \code{3x3} symmetric positive definite matrix.
#' @param isigma2 A \code{3x3} symmetric positive definite matrix.
#' @param neutral_variance The neutral element of the Bayes space (default:
#'   1e5).
#'
#' @return A positive scalar measuring the distance between the two input
#'   Gaussians in Bayes geometry.
#' @export
#'
#' @examples
#' D1 <- diag(c(1.7, 0.3, 0.1)) * 1e-3
#' D2 <- diag(c(0.3, 1.7, 0.1)) * 1e-3
#' bayes_distance(D1, D2)
bayes_distance <- function(isigma1, isigma2, neutral_variance = neutral_element) {
  isigma1 <- isigma1 * neutral_variance
  isigma2 <- isigma2 * neutral_variance
  diff <- isigma1 - isigma2
  sum(diff^2) / 2
}
