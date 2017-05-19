#' Addition in Bayes space
#'
#' @param isigma1 A \code{3x3} symmetric positive definite matrix.
#' @param isigma2 A \code{3x3} symmetric positive definite matrix.
#' @param neutral The neutral element of the Bayes space (default:
#'   \code{diag(3e-3, 3L)}).
#'
#' @return A \code{3x3} symmetric positive definite matrix.
#' @export
#'
#' @examples
#' D1 <- diag(c(1.7, 0.3, 0.1)) * 1e-3
#' D2 <- diag(c(0.3, 1.7, 0.1)) * 1e-3
#' bayes_add(D1, D2)
bayes_add <- function(isigma1, isigma2, neutral = neutral_element) {
  isigma1 + isigma2 - neutral
}

#' Scalar multiplication in Bayes space
#'
#' @param isigma A \code{3x3} symmetric positive definite matrix.
#' @param alpha A scalar.
#' @param neutral The neutral element of the Bayes space (default:
#'   \code{diag(3e-3, 3L)}).
#'
#' @return A \code{3x3} symmetric positive definite matrix.
#' @export
#'
#' @examples
#' D <- diag(c(1.7, 0.3, 0.1)) * 1e-3
#' bayes_multiply(D, 2)
#' bayes_add(D, D)
bayes_multiply <- function(isigma, alpha, neutral = neutral_element) {
  alpha * isigma - (alpha - 1) * neutral
}

#' Mean element in Bayes space
#'
#' @param isigmaList A list of \code{3x3} symmetric positive definite matrices.
#' @param weights A numeric vector specifying the weight of each observation
#'   (default: equiprobability).
#' @param lambda A scalar between 0 and 1 specifying how much trust do we have
#'   in the reference measure, which is assigned a weight of \code{1 - lambda}
#'   (default: 1).
#' @param neutral The neutral element of the Bayes space (default:
#'   \code{diag(3e-3, 3L)}.
#'
#' @return A \code{3x3} symmetric positive definite matrix.
#' @export
#'
#' @examples
#' D <- list(D1 = diag(c(1.7, 0.3, 0.1)) * 1e-3,
#'           D2 = diag(c(0.3, 1.7, 0.1)) * 1e-3)
#' bayes_mean(D)
bayes_mean <- function(isigmaList,
                       weights = rep(1 / length(isigmaList), length(isigmaList)),
                       lambda = 1,
                       neutral = neutral_element) {
  weights <- weights / sum(weights)
  isigmaList %>%
    purrr::map2(weights * lambda, bayes_multiply, neutral = neutral) %>%
    purrr::reduce(bayes_add, neutral = neutral, .init = neutral)
}

#' Distance in Bayes space
#'
#' @param isigma1 A \code{3x3} symmetric positive definite matrix.
#' @param isigma2 A \code{3x3} symmetric positive definite matrix.
#' @param neutral The neutral element of the Bayes space (default:
#'   \code{diag(3e-3, 3L)}).
#'
#' @return A positive scalar measuring the distance between the two input
#'   Gaussians in Bayes geometry.
#' @export
#'
#' @examples
#' D1 <- diag(c(1.7, 0.3, 0.1)) * 1e-3
#' D2 <- diag(c(0.3, 1.7, 0.1)) * 1e-3
#' bayes_distance(D1, D2)
bayes_distance <- function(isigma1, isigma2, neutral = neutral_element) {
  # Needs recalculation to accomodate a non-identity neutral element
  isigma1 <- isigma1 * neutral
  isigma2 <- isigma2 * neutral
  diff <- isigma1 - isigma2
  sum(diff^2) / 2
}
