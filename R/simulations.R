#' Robustness of log-Euclidean and Bayes geometries to MRI-induced noise
#'
#' @param tensor A reference tensor.
#' @param n Sample size (default: \code{8L}).
#' @param B Number of independent noisy samples (default: \code{1000L}).
#' @param N Number of healthy subjects used to build hypothetical template
#'   (default: \code{20L}).
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
robustness_analysis <- function(tensor, n = 8L, B = 1000L, N = 20L, seed = NULL) {
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
    multidplyr::cluster_copy(vec2mat) %>%
    multidplyr::cluster_copy(mat2vec) %>%
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
        weights = weights
      ))
    ) %>%
    dplyr::collect() %>%
    dplyr::ungroup() %>%
    tidyr::unnest()
  list(data = data, weights = weights)
}

single_run <- function(tensor, N, rho, n, B, weights = rep(1 / n, n)) {
  tibble::tibble(Replicate = paste0("Rep", seq_len(B))) %>%
    dplyr::mutate(
      Estimates = purrr::map(Replicate, ~ average_estimators(
        tensor = tensor,
        N = N,
        rho = rho,
        n = n,
        weights = weights
      )),
      Euclidean_Estimate = purrr::map(Estimates, "Euclidean"),
      LogEuclidean_Estimate = purrr::map(Estimates, "LogEuclidean"),
      Bayes10_Estimate = purrr::map(Estimates, "Bayes10"),
      Bayes20_Estimate = purrr::map(Estimates, "Bayes20"),
      Bayes30_Estimate = purrr::map(Estimates, "Bayes30"),
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
      Bayes10_Microstructure = purrr::map(Bayes10_Estimate,
                                               microstructure_distance, tensor),
      Bayes10_Diffusivity = purrr::map_dbl(Bayes10_Microstructure, "diffusivity"),
      Bayes10_Radius = purrr::map_dbl(Bayes10_Microstructure, "radius"),
      Bayes10_Direction = purrr::map_dbl(Bayes10_Microstructure, "direction"),
      Bayes20_Microstructure = purrr::map(Bayes20_Estimate,
                                        microstructure_distance, tensor),
      Bayes20_Diffusivity = purrr::map_dbl(Bayes20_Microstructure, "diffusivity"),
      Bayes20_Radius = purrr::map_dbl(Bayes20_Microstructure, "radius"),
      Bayes20_Direction = purrr::map_dbl(Bayes20_Microstructure, "direction"),
      Bayes30_Microstructure = purrr::map(Bayes30_Estimate,
                                        microstructure_distance, tensor),
      Bayes30_Diffusivity = purrr::map_dbl(Bayes30_Microstructure, "diffusivity"),
      Bayes30_Radius = purrr::map_dbl(Bayes30_Microstructure, "radius"),
      Bayes30_Direction = purrr::map_dbl(Bayes30_Microstructure, "direction")
    ) %>%
    dplyr::select(9:11, 13:15, 17:19, 21:23, 25:27) %>%
    dplyr::summarise_all(mean) %>%
    tidyr::gather(Tmp, MSE) %>%
    tidyr::separate(Tmp, c("Space", "Metric"))
}

average_estimators <- function(tensor, N, rho, n, weights = rep(1 / n, n)) {
  tensors <- rmrice(tensor, n = n, rho = rho)
  tensor_ref <- rmrice(tensor, n = 1, rho = rho / sqrt(N))[[1]]
  list(
    Euclidean = mean_euclidean(tensors, weights),
    LogEuclidean = mean_log_euclidean(tensors, weights),
    Bayes10 = mean_bayes(tensors, weights, shrinkage = exp(-(1 / rho) / 10), tensor_ref = tensor_ref),
    Bayes20 = mean_bayes(tensors, weights, shrinkage = exp(-(1 / rho) / 20), tensor_ref = tensor_ref),
    Bayes30 = mean_bayes(tensors, weights, shrinkage = exp(-(1 / rho) / 30), tensor_ref = tensor_ref)
  )
}

#' @export
estimation_analysis <- function(tensor, B = 1000L, N = 20L, seed = NULL) {
  set.seed(seed)
  cl <- multidplyr::create_cluster(cores = 4L) %>%
    multidplyr::cluster_library("tidyverse") %>%
    multidplyr::cluster_copy(single_run_sim2) %>%
    multidplyr::cluster_copy(rmrice) %>%
    multidplyr::cluster_copy(B) %>%
    multidplyr::cluster_copy(N) %>%
    multidplyr::cluster_copy(tensor) %>%
    multidplyr::cluster_copy(vec2mat) %>%
    multidplyr::cluster_copy(mat2vec) %>%
    multidplyr::cluster_copy(exp_tensor) %>%
    multidplyr::cluster_copy(log_tensor) %>%
    multidplyr::cluster_copy(dist_euclidean)
  tibble::tibble(Sigma = seq(2, 16, by = 2) / 100) %>%
    multidplyr::partition(cluster = cl) %>%
    dplyr::mutate(
      data = purrr::map(Sigma, ~ single_run_sim2(
        tensor = tensor,
        N = N,
        rho = .x,
        B = B
      ))
    ) %>%
    dplyr::collect() %>%
    dplyr::ungroup() %>%
    tidyr::unnest()
}

single_run_sim2 <- function(tensor, N, rho, B) {
  signals <- S0 * exp(-bval * diag(bvecs %*% tensor %*% t(bvecs)))
  tibble::tibble(Replicate = paste0("Rep", seq_len(B))) %>%
    dplyr::mutate(
      Signals = signals %>% purrr::map(rrice, n = B, sigma = rho * S0) %>%
        purrr::transpose() %>%
        purrr::simplify_all(),

      GaussianSignals_Estimate = purrr::map(
        Signals,
        estimate_tensor,
        tensor_init = tensor,
        input_type = "original",
        sigma = NULL
      ),
      RicianSignals_Estimate = purrr::map(
        Signals,
        estimate_tensor,
        tensor_init = tensor,
        input_type = "original",
        sigma = rho * S0
      ),
      LogSignalsOLS_Estimate = purrr::map(
        Signals,
        estimate_tensor,
        tensor_init = tensor,
        input_type = "log",
        wls = FALSE
      ),
      LogSignalsWLS_Estimate = purrr::map(
        Signals,
        estimate_tensor,
        tensor_init = tensor,
        input_type = "log",
        wls = TRUE
      ),

      GaussianSignals_Microstructure = purrr::map(
        GaussianSignals_Estimate,
        microstructure_distance,
        tensor
      ),
      GaussianSignals_Diffusivity = purrr::map_dbl(GaussianSignals_Microstructure, "diffusivity"),
      GaussianSignals_Radius = purrr::map_dbl(GaussianSignals_Microstructure, "radius"),
      GaussianSignals_Direction = purrr::map_dbl(GaussianSignals_Microstructure, "direction"),

      RicianSignals_Microstructure = purrr::map(
        RicianSignals_Estimate,
        microstructure_distance,
        tensor
      ),
      RicianSignals_Diffusivity = purrr::map_dbl(RicianSignals_Microstructure, "diffusivity"),
      RicianSignals_Radius = purrr::map_dbl(RicianSignals_Microstructure, "radius"),
      RicianSignals_Direction = purrr::map_dbl(RicianSignals_Microstructure, "direction"),

      LogSignalsOLS_Microstructure = purrr::map(
        LogSignalsOLS_Estimate,
        microstructure_distance,
        tensor
      ),
      LogSignalsOLS_Diffusivity = purrr::map_dbl(LogSignalsOLS_Microstructure, "diffusivity"),
      LogSignalsOLS_Radius = purrr::map_dbl(LogSignalsOLS_Microstructure, "radius"),
      LogSignalsOLS_Direction = purrr::map_dbl(LogSignalsOLS_Microstructure, "direction"),

      LogSignalsWLS_Microstructure = purrr::map(
        LogSignalsWLS_Estimate,
        microstructure_distance,
        tensor
      ),
      LogSignalsWLS_Diffusivity = purrr::map_dbl(LogSignalsWLS_Microstructure, "diffusivity"),
      LogSignalsWLS_Radius = purrr::map_dbl(LogSignalsWLS_Microstructure, "radius"),
      LogSignalsWLS_Direction = purrr::map_dbl(LogSignalsWLS_Microstructure, "direction")
    ) %>%
    dplyr::select(8:10, 12:14, 16:18, 20:22) %>%
    dplyr::summarise_all(mean) %>%
    tidyr::gather(Tmp, MSE) %>%
    tidyr::separate(Tmp, c("Space", "Metric"))
}
