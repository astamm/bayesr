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
#' mat2vec(D)
mat2vec <- function(tensor, twice = FALSE) {
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
#' vec2mat(v)
vec2mat <- function(vector, twice = FALSE) {
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
    purrr::map(estimate_tensor, tensor_init = tensor)
}

is_pd <- function(tensor) {
  eigen(tensor, TRUE, TRUE)$values[3] > sqrt(.Machine$double.eps)
}

estimate_tensor <- function(signals, tensor_init, input_type = "original", wls = FALSE, sigma = NULL) {
  if (input_type == "original") {
    # Estimation based on original signals
    if (is.null(sigma)) { # Gaussian Likelihood Maximization
      opt <- nloptr::newuoa(
        x0 = c(log(S0), mat2vec(log_tensor(tensor_init))),
        fn = likelihood_gaussian,
        observed_signals = signals,
        nl.info = TRUE,
        control = list(xtol_rel = 1e-8, maxeval = 1e5)
      )
      return(exp_tensor(vec2mat(opt$par[-1])))
    } else { # Rician Likelihood Maximization
      opt <- nloptr::newuoa(
        x0 = c(log(S0), log(sigma^2), mat2vec(log_tensor(tensor_init))),
        fn = likelihood_rician,
        observed_signals = signals,
        nl.info = TRUE,
        control = list(xtol_rel = 1e-8, maxeval = 1e5)
      )
      return(exp_tensor(vec2mat(opt$par[-(1:2)])))
    }
  }

  # Estimation based on log-signals
  log_signals <- log(signals)

  if (wls) { # Weighted Least Squares
    logS0 <- 0
    logTmp <- 1
    while (abs(logS0 - logTmp) > 1e-3) {
      logS0 <- logTmp
      weights <- exp(as.numeric(design_matrix %*% estimation_matrix %*% log_signals))
      W2 <- diag(weights^2)
      iM <- solve(t(design_matrix) %*% W2 %*% design_matrix)
      tensor <- iM %*% t(design_matrix) %*% W2 %*% log_signals
      logTmp <- tensor[1]
      log_signals <- as.numeric(design_matrix %*% tensor)
    }
  } else # Original Least Squares
    tensor <- estimation_matrix %*% log_signals

  vec2mat(tensor[-1])
}

likelihood_gaussian <- function(x, observed_signals) {
  # Parameters are log(S0) [1] and log(D) [2:7]
  S0 <- exp(x[1])
  tensor <- exp_tensor(vec2mat(x[-1]))
  predicted_signals <- S0 * exp(-bval * diag(bvecs %*% tensor %*% t(bvecs)))
  sum((observed_signals - predicted_signals)^2)
}

likelihood_rician <- function(x, observed_signals) {
  # Parameters are log(S0) [1], log(\sigma^2) [2] and log(D) [3:8]
  S0 <- exp(x[1])
  sigma2 <- exp(x[2])
  tensor <- exp_tensor(vec2mat(x[-(1:2)]))
  predicted_signals <- S0 * exp(-bval * diag(bvecs %*% tensor %*% t(bvecs)))
  -2 * sum(
    log(observed_signals) - log(sigma2) -
      (observed_signals - predicted_signals)^2 / (2 * sigma2) +
      log(besselI(observed_signals * predicted_signals / sigma2, nu = 0, expon.scaled = TRUE))
  )
}

#' Distances on brain microstructure
#'
#' @param tensor1 A symmetric definite positive 3x3 matrix.
#' @param tensor2 A symmetric definite positive 3x3 matrix.
#'
#' @return A list of distances for microstructure parameters.
#' @export
#'
#' @examples
#' D1 <- diag(c(1.7, 0.3, 0.1)) * 1e-3
#' D2 <- diag(c(0.3, 1.7, 0.1)) * 1e-3
#' microstructure_distance(D1, D2)
microstructure_distance <- function(tensor1, tensor2) {
  eig1 <- eigen(tensor1, symmetric = TRUE)
  eig2 <- eigen(tensor2, symmetric = TRUE)
  l1 <- eig1$values[1]
  l2 <- eig2$values[1]
  r1 <- mean(eig1$values[2:3])
  r2 <- mean(eig2$values[2:3])
  dir1 <- eig1$vectors[, 1]
  dir2 <- eig2$vectors[, 1]
  list(
    diffusivity = (l1 - l2)^2,
    radius = (r1 - r2)^2,
    direction = acos(abs(sum(dir1 * dir2)))^2
  )
}

log_tensor <- function(x) {
  eig <- eigen(x, symmetric = TRUE)
  val <- log(eig$values)
  eig$vectors %*% diag(val) %*% t(eig$vectors)
}

exp_tensor <- function(x) {
  eig <- eigen(x, symmetric = TRUE)
  val <- exp(eig$values)
  eig$vectors %*% diag(val) %*% t(eig$vectors)
}

uniform_weights <- function(tensor_list) {
  len <- length(tensor_list)
  rep(1 / len, len)
}

#' Determinant of a symmetric matrix
#'
#' @param x A symmetric matrix.
#'
#' @return The determinant of the input symmetric matrix.
#' @export
#'
#' @examples
#' M <- matrix(rnorm(9), 3, 3)
#' det_matrix(M)
det_matrix <- function(x) {
  prod(eigen(x, TRUE, TRUE)$values)
}
