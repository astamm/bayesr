#' Robustness of log-Euclidean and Bayes geometries to MRI-induced noise
#'
#' @param tensor A reference tensor.
#' @param n Sample size (default: \code{8L}).
#' @param B Number of independent noisy samples (default: \code{1000L}).
#'
#' @return A \code{\link[tibble]{tibble}} with simulation results.
#' @export
#'
#' @examples
#' arrow <- ggplot2::arrow(length = ggplot2::unit(0.4, "cm"), type = "closed")
#'
#' refIsotropicTensor <- diag(3e-3, 3L)
#' data_isotropic <- robustness_analysis(refIsotropicTensor, B = 100L)
#'
#' data_isotropic %>%
#'  ggplot2::ggplot(ggplot2::aes(x = Sigma, y = MSE, col = Space)) +
#'  ggplot2::geom_point() +
#'  ggplot2::geom_line() +
#'  ggplot2::theme_minimal() +
#'  ggplot2::theme(legend.position = "top",
#'                 axis.line = ggplot2::element_line(arrow = arrow)) +
#'  ggplot2::facet_grid(Metric ~ ., scales = "free") +
#'  ggplot2::scale_x_continuous(labels = scales::percent) +
#'  ggplot2::scale_y_log10()
#'
#' refAnisotropicTensor <- diag(c(1.71e-3, 3e-4, 1e-4))
#' data_fascicles <- tibble::tibble()
#' theta <- pi * c(0, 1/6, 1/4, 1/3, 1/2)
#' for (a in theta) {
#'  R <- rbind(
#'    c(cos(a), sin(a), 0),
#'    c(-sin(a), cos(a), 0),
#'    c(0, 0, 1)
#'  )
#'  ref_tmp <- R %*% refAnisotropicTensor %*% t(R)
#'  data_tmp <- robustness_analysis(ref_tmp, B = 100L)
#'  data_fascicles <- dplyr::bind_rows(
#'    data_fascicles,
#'    data_tmp %>% dplyr::mutate(Angle = round(a, 4L))
#'  )
#' }
#'
#' data_fascicles %>%
#'   ggplot2::ggplot(ggplot2::aes(x = Sigma, y = MSE, col = Space)) +
#'   ggplot2::geom_point() +
#'   ggplot2::geom_line() +
#'   ggplot2::theme_minimal() +
#'   ggplot2::theme(legend.position = "top",
#'                  axis.line = ggplot2::element_line(arrow = arrow)) +
#'   ggplot2::facet_grid(Metric ~ Angle, scales = "free") +
#'   ggplot2::scale_x_continuous(labels = scales::percent) +
#'   ggplot2::scale_y_log10()
robustness_analysis <- function(tensor, n = 8L, B = 1000L) {
  tibble::tibble(Sigma = seq(2, 16, by = 2) / 100) %>%
    dplyr::mutate(
      data = purrr::map(Sigma, single_run, tensor = tensor, n = n, B = B)
    ) %>%
    tidyr::unnest()
}

single_run <- function(tensor, rho, n, B) {
  tibble::tibble(Replicate = paste0("Rep", seq_len(B))) %>%
    dplyr::mutate(
      Estimates = purrr::map(Replicate, average_estimators, tensor = tensor,
                             rho = rho, n = n),
      LogEuclidean_Estimate = purrr::map(Estimates, "LogEuclidean"),
      Bayes_Estimate = purrr::map(Estimates, "Bayes"),
      LogEuclidean_LogEuclidean = purrr::map_dbl(LogEuclidean_Estimate,
                                                 log_euclidean_distance, tensor),
      LogEuclidean_Bayes = purrr::map_dbl(LogEuclidean_Estimate, bayes_distance,
                                          tensor),
      Bayes_LogEuclidean = purrr::map_dbl(Bayes_Estimate, log_euclidean_distance,
                                          tensor),
      Bayes_Bayes = purrr::map_dbl(Bayes_Estimate, bayes_distance, tensor)
    ) %>%
    dplyr::select(5:8) %>%
    dplyr::summarise_all(mean) %>%
    tidyr::gather(Tmp, MSE) %>%
    tidyr::separate(Tmp, c("Space", "Metric"))
}

average_estimators <- function(tensor, rho, n, id = "Rep1") {
  tensors <- rmrice(tensor, n = n, rho = rho)
  list(LogEuclidean = log_euclidean_mean(tensors), Bayes = bayes_mean(tensors))
}
