mean_euclidean <- function(tensor_list, weights = uniform_weights(tensor_list)) {
  weights <- weights / sum(weights)
  tensor_list %>%
    purrr::map2(weights, `*`) %>%
    purrr::reduce(`+`)
}

dist_euclidean <- function(tensor1, tensor2) {
  norm(tensor1 - tensor2, "F")^2
}
