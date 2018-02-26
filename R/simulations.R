#' Robustness of log-Euclidean and Bayes geometries to MRI-induced noise
#'
#' @param tensor A reference tensor.
#' @param n Sample size (default: \code{8L}).
#' @param B Number of independent noisy samples (default: \code{1000L}).
#' @param N Number of healthy subjects used to build hypothetical template
#'   (default: \code{20L}).
#' @param beta Scaling parameter for defining SNR weight for reference measure
#'   (default: \code{20}).
#' @param seed Seed for random number generation (default: clock time).
#'
#' @return A \code{\link[tibble]{tibble}} with simulation results.
#' @export
#'
#' @examples
#' arrow <- ggplot2::arrow(length = ggplot2::unit(0.4, "cm"), type = "closed")
#'
#' refIsotropicTensor <- diag(3e-3, 3L)
#' data_isotropic <- robustness_analysis(refIsotropicTensor, B = 100L, seed = 1234)
#'
#' data_isotropic$data %>%
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
#'  tmp <- robustness_analysis(ref_tmp, B = 100L, seed = 1234)
#'  data_fascicles <- dplyr::bind_rows(
#'    data_fascicles,
#'    tmp$data %>% dplyr::mutate(Angle = round(a, 4L))
#'  )
#'}
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
robustness_analysis <- function(tensor, n = 8L, B = 1000L, N = 20L, beta = 20, seed = NULL) {
  set.seed(seed)
  weights <- runif(n)
  cl <- multidplyr::create_cluster(cores = 4L) %>%
    multidplyr::cluster_library("tidyverse") %>%
    multidplyr::cluster_copy(single_run) %>%
    multidplyr::cluster_copy(average_estimators) %>%
    multidplyr::cluster_copy(rmrice) %>%
    multidplyr::cluster_copy(B) %>%
    multidplyr::cluster_copy(n) %>%
    multidplyr::cluster_copy(N) %>%
    multidplyr::cluster_copy(tensor) %>%
    multidplyr::cluster_copy(weights) %>%
    multidplyr::cluster_copy(beta) %>%
    multidplyr::cluster_copy(as_tensor) %>%
    multidplyr::cluster_copy(as_vector) %>%
    multidplyr::cluster_copy(exp_tensor) %>%
    multidplyr::cluster_copy(log_tensor)
  data <- tibble::tibble(Sigma = seq(2, 16, by = 2) / 100) %>%
    multidplyr::partition(cluster = cl) %>%
    dplyr::mutate(
      data = purrr::map(Sigma, ~ single_run(
        tensor = tensor,
        N = N,
        rho = .x,
        n = n,
        B = B,
        weights = weights,
        shrinkage = exp(-(1 / .x) / beta)
      ))
    ) %>%
    dplyr::collect() %>%
    dplyr::ungroup() %>%
    tidyr::unnest()
  list(data = data, weights = weights)
}

single_run <- function(tensor, N, rho, n, B, weights = rep(1 / n, n), shrinkage = 0) {
  tibble::tibble(Replicate = paste0("Rep", seq_len(B))) %>%
    dplyr::mutate(
      Estimates = purrr::map(Replicate, ~ average_estimators(
        tensor = tensor,
        N = N,
        rho = rho,
        n = n,
        weights = weights,
        shrinkage = shrinkage
      )),
      Euclidean_Estimate = purrr::map(Estimates, "Euclidean"),
      LogEuclidean_Estimate = purrr::map(Estimates, "LogEuclidean"),
      Bayes_Estimate = purrr::map(Estimates, "Bayes"),
      Euclidean_Microstructure = purrr::map(Euclidean_Estimate,
                                            microstructure_distance, tensor),
      Euclidean_Diffusivity = purrr::map_dbl(Euclidean_Microstructure, "diffusivity"),
      Euclidean_Radius = purrr::map_dbl(Euclidean_Microstructure, "radius"),
      Euclidean_Direction = purrr::map_dbl(Euclidean_Microstructure, "direction"),
      LogEuclidean_Microstructure = purrr::map(LogEuclidean_Estimate,
                                               microstructure_distance, tensor),
      LogEuclidean_Diffusivity = purrr::map_dbl(LogEuclidean_Microstructure, "diffusivity"),
      LogEuclidean_Radius = purrr::map_dbl(LogEuclidean_Microstructure, "radius"),
      LogEuclidean_Direction = purrr::map_dbl(LogEuclidean_Microstructure, "direction"),
      Bayes_Microstructure = purrr::map(Bayes_Estimate,
                                               microstructure_distance, tensor),
      Bayes_Diffusivity = purrr::map_dbl(Bayes_Microstructure, "diffusivity"),
      Bayes_Radius = purrr::map_dbl(Bayes_Microstructure, "radius"),
      Bayes_Direction = purrr::map_dbl(Bayes_Microstructure, "direction")
    ) %>%
    dplyr::select(7:9, 11:13, 15:17) %>%
    dplyr::summarise_all(mean) %>%
    tidyr::gather(Tmp, MSE) %>%
    tidyr::separate(Tmp, c("Space", "Metric"))
}

average_estimators <- function(tensor, N, rho, n, weights = rep(1 / n, n), shrinkage = 0) {
  tensors <- rmrice(tensor, n = n, rho = rho)
  tensor_ref <- rmrice(tensor, n = 1, rho = rho / sqrt(N))[[1]]
  list(
    Euclidean = mean_euclidean(tensors, weights),
    LogEuclidean = mean_log_euclidean(tensors, weights),
    Bayes = mean_bayes(tensors, weights, shrinkage = shrinkage, tensor_ref = tensor_ref)
  )
}
