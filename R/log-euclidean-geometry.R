mean_log_euclidean <- function(tensor_list, weights = uniform_weights(tensor_list)) {
  weights <- weights / sum(weights)
  tensor_list %>%
    purrr::map(log_tensor) %>%
    purrr::map2(weights, `*`) %>%
    purrr::reduce(`+`) %>%
    exp_tensor()
}

dist_log_euclidean <- function(tensor1, tensor2) {
  norm(log_tensor(tensor1) - log_tensor(tensor2), type = "F")^2
}
