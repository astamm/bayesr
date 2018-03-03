#' Addition in Bayes space
#'
#' @param tensor1 A \code{3x3} symmetric matrix.
#' @param tensor2 A \code{3x3} symmetric matrix.
#' @param tensor_ref The SDP tensor associated to the reference measure of the
#'   Bayes space.
#' @param do_projection A boolean specifying whether the result of the operation
#'   should be projected onto the space of diffusion signals (default not to).
#'
#' @return A \code{3x3} symmetric matrix.
#' @export
#'
#' @examples
#' D1 <- diag(c(1.7, 0.3, 0.1)) * 1e-3
#' D2 <- diag(c(0.3, 1.7, 0.1)) * 1e-3
#' Dr <- diag(c(1.7, 0.3, 0.1)) * 1e-3
#' add_bayes(D1, D2, Dr)
add_bayes <- function(tensor1, tensor2, tensor_ref, do_projection = FALSE) {
  res <- tensor1 + tensor2 - tensor_ref

  if (do_projection & !is_pd(res))
    return(project(res, tensor_ref))

  res
}

#' Scalar multiplication in Bayes space
#'
#' @param tensor A \code{3x3} symmetric matrix.
#' @param alpha A scalar.
#' @param tensor_ref The SDP tensor associated to the reference measure of the
#'   Bayes space.
#' @param do_projection A boolean specifying whether the result of the operation
#'   should be projected onto the space of diffusion signals (default not to).
#'
#' @return A \code{3x3} symmetric matrix.
#' @export
#'
#' @examples
#' D <- diag(c(1.7, 0.3, 0.1)) * 1e-3
#' multiply_bayes(D, 2)
multiply_bayes <- function(tensor, alpha, tensor_ref, do_projection = FALSE) {
  res <- alpha * tensor + (1 - alpha) * tensor_ref

  if (do_projection & !is_pd(res))
    return(project(res, tensor_ref))

  res
}

#' Mean element in Bayes space
#'
#' @param tensor_list A list of \code{3x3} symmetric matrices.
#' @param weights A numeric vector specifying the weight of each observation
#'   (default: equiprobability).
#' @param shrinkage A scalar between 0 and 1 specifying how much do we want to
#'   shrink the arithmetic mean towards the reference tensor (default: 0).
#' @param tensor_ref The SDP tensor associated to the reference measure of the
#'   Bayes space.
#' @param do_projection A boolean specifying whether the result of the operation
#'   should be projected onto the space of diffusion signals (default not to).
#'
#' @return A \code{3x3} symmetric matrix.
#' @export
#'
#' @examples
#' D <- list(D1 = diag(c(1.7, 0.3, 0.1)) * 1e-3,
#'           D2 = diag(c(0.3, 1.7, 0.1)) * 1e-3)
#' mean_bayes(D)
mean_bayes <- function(tensor_list,
                       weights = uniform_weights(tensor_list),
                       shrinkage = 0,
                       tensor_ref,
                       do_projection = FALSE) {
  weights <- weights / sum(weights) * (1 - shrinkage)
  res <- tensor_list %>%
    purrr::map2(weights, multiply_bayes, tensor_ref = tensor_ref, do_projection = do_projection) %>%
    purrr::reduce(add_bayes, tensor_ref = tensor_ref, do_projection = do_projection, .init = tensor_ref)

  if (do_projection & !is_pd(res))
    return(project(res, tensor_ref))

  res
}

#' Distance in Bayes space
#'
#' @param tensor1 A \code{3x3} symmetric matrix.
#' @param tensor2 A \code{3x3} symmetric matrix.
#' @param iR The inverse of the SDP tensor associated to the reference measure
#'   of the Bayes space.
#'
#' @return A positive scalar.
#' @export
#'
#' @examples
#' D1 <- diag(c(1.7, 0.3, 0.1)) * 1e-3
#' D2 <- diag(c(0.3, 1.7, 0.1)) * 1e-3
#' dist_bayes(D1, D2)
dist_bayes <- function(tensor1, tensor2, iR) {
  M <- (tensor1 - tensor2) %*% iR
  vals <- eigen(M, only.values = TRUE)$values
  # sum(vals^2) #* sqrt(pi^3 * det_matrix(iR))
  sum(M * t(M)) * sqrt(pi^3 * det_matrix(iR))
}

cost <- function(x, S, iR, distance = "bayes") {
  DT <- vec2mat(x)
  DT <- exp_tensor(DT)
  if (distance != "bayes")
    return(dist_euclidean(S, DT) + dist_euclidean(solve(iR), DT))
  dist_bayes(S, DT, iR) + dist_bayes(solve(iR), DT, iR)
}

#' Orthogonal projection onto diffusion signal space
#'
#' @param x A symmetric 3x3 matrix resulting from a series of operations in
#'   Bayes space.
#' @param reference_tensor The SDP tensor defining the reference measure of the
#'   Bayes space.
#' @param distance A string indicating whether the projection should be done in
#'   Bayes geometry (default) or in Euclidean geometry.
#'
#' @return An SDP tensor.
#' @export
#'
#' @examples
#' R <- diag(c(1.7, 0.3, 0.1), 3)
#' Rn <- matrix(mvtnorm::rmvnorm(1, c(R)), 3, 3)
#' det_matrix(Rn)
#' Rfro <- project(Rn, R, distance = "frobenius")
#' Rbay <- project(Rn, R, distance = "bayes")
#' eigen(Rn, TRUE, TRUE)$values
#' eigen(Rfro, TRUE, TRUE)$values
#' eigen(Rbay, TRUE, TRUE)$values
project <- function(x, reference_tensor, distance = "bayes") {
  logR <- log_tensor(reference_tensor)
  invR <- solve(reference_tensor)
  opt <- nloptr::newuoa(
    x0 = mat2vec(logR),
    fn = cost,
    S = x, iR = invR, distance = distance,
    nl.info = TRUE,
    control = list(xtol_rel = 1e-8, maxeval = 1e5)
  )
  if (opt$value < 0)
    print(paste(distance, opt$value))
  op <- opt$par
  logT <- vec2mat(op)
  exp_tensor(logT)
}
