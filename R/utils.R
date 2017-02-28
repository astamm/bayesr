#' @importFrom dplyr %>%
#' @export
dplyr::`%>%`

rrice <- function (n, mu, sigma)
{
  theta <- 1
  x <- rnorm(n, mean = mu * cos(theta), sd = sigma)
  y <- rnorm(n, mean = mu * sin(theta), sd = sigma)
  sqrt(x^2 + y^2)
}

#' Tensor-To-Vector Representation
#'
#' @param tensor A 3x3 symmetric postive definite matrix.
#' @param twice A boolean that says whether off-diagonal tensor elements were
#'   doubled in vector representation (default: \code{FALSE}).
#'
#' @return A vector of size 6 storing the 6 unique components of the input tensor.
#' @export
#'
#' @examples
#' D <- cbind(c(1, 2, 4), c(2, 3, 5), c(4, 5, 6))
#' as_vector(D)
as_vector <- function(tensor, twice = FALSE) {
  stopifnot(is.matrix(tensor))
  if (!is.matrix(tensor))
    stop("Input should be a matrix")
  else {
    if (any(dim(tensor) != c(3L, 3L)))
      stop("Input should be of dimension 3x3")
  }

  if (twice) {
    tensor[2, 1] <- 2 * tensor[2, 1]
    tensor[3, 1] <- 2 * tensor[3, 1]
    tensor[3, 2] <- 2 * tensor[3, 2]
  }

  c(tensor[1, 1], tensor[2, 1], tensor[2, 2],
    tensor[3, 1], tensor[3, 2], tensor[3, 3])
}

#' Vector-To-Tensor Representation
#'
#' @param vector A vector of size 6 storing the 6 unique components of a tensor.
#' @param twice A boolean that says whether off-diagonal tensor elements were
#'   doubled in vector representation (default: \code{FALSE}).
#'
#' @return A 3x3 symmetric postive definite matrix as output tensor.
#' @export
#'
#' @examples
#' v <- seq_len(6L)
#' as_tensor(v)
as_tensor <- function(vector, twice = FALSE) {
  if (!is.numeric(vector))
    stop("Input should be a numeric vector")
  else {
    if (length(vector) != 6L)
      stop("Input should be of dimension 6")
  }

  xx <- vector[1]
  yx <- vector[2]
  yy <- vector[3]
  zx <- vector[4]
  zy <- vector[5]
  zz <- vector[6]

  if (twice) {
    yx <- yx / 2
    zx <- zx / 2
    zy <- zy / 2
  }

  cbind(c(xx, yx, zx), c(yx, yy, zy), c(zx, zy, zz))
}

#' MRI Noise Generator
#'
#' @param tensor A symmetric definite positive 3x3 matrix.
#' @param n Number of noisy independent versions of the input tensor to be
#'   generated (default: 1).
#' @param rho Relative standard deviation of the MRI noise (default: 0.1).
#'
#' @return A list of \code{n} noisy independent
#'   versions of the input tensor.
#' @export
#'
#' @examples
#' D <- diag(c(1.7, 0.3, 0.1)) * 1e-3
#' rmrice(D, 5L, 0.05)
rmrice <- function(tensor, n = 1L, rho = 0.1) {
  signals <- S0 * exp(-bval * diag(bvecs %*% tensor %*% t(bvecs)))
  signals %>% purrr::map(rrice, n = n, sigma = rho * S0) %>%
    purrr::transpose() %>%
    purrr::simplify_all() %>%
    purrr::map(estimate_tensor)
}

estimate_tensor <- function(signals) {
  tensor <- (estimation_matrix %*% log(signals))[2:7]
  tensor <- as_tensor(tensor)
  as.matrix(Matrix::nearPD(tensor)$mat)
}
