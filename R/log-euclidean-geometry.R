log_euclidean_mean <- function(sigmaList, weights = rep(1 / length(sigmaList), length(sigmaList))) {
  weights <- weights / sum(weights)
  sigmaList %>%
    purrr::map(log_tensor) %>%
    purrr::map2(weights, `*`) %>%
    purrr::reduce(`+`) %>%
    exp_tensor()
}

log_euclidean_distance <- function(sigma1, sigma2) {
  diff <- log_tensor(sigma1) - log_tensor(sigma2)
  sum(diff^2)
}
