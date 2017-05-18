euclidean_mean <- function(sigmaList, weights = rep(1 / length(sigmaList), length(sigmaList))) {
  sigmaList %>%
    purrr::map2(weights, `*`) %>%
    purrr::reduce(`+`)
}

euclidean_distance <- function(sigma1, sigma2) {
  norm(sigma1 - sigma2, "F")^2
}
