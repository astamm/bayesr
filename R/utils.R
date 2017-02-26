#' Tensor-To-Vector Representation
#'
#' @param tensor A 3x3 symmetric postive definite matrix.
#'
#' @return A vector of size 6 storing the 6 unique components of the input tensor.
#' @export
#'
#' @examples
#' D <- cbind(c(1, 2, 4), c(2, 3, 5), c(4, 5, 6))
#' as_vector(D)
as_vector <- function(tensor) {
  stopifnot(is.matrix(tensor))
  if (!is.matrix(tensor))
    stop("Input should be a matrix")
  else {
    if (any(dim(tensor) != c(3L, 3L)))
      stop("Input should be of dimension 3x3")
  }

  c(tensor[1, 1], tensor[2, 1], tensor[2, 2],
    tensor[3, 1], tensor[3, 2], tensor[3, 3])
}

#' Vector-To-Tensor Representation
#'
#' @param vector A vector of size 6 storing the 6 unique components of a tensor.
#'
#' @return A 3x3 symmetric postive definite matrix as output tensor.
#' @export
#'
#' @examples
#' v <- seq_len(6L)
#' as_tensor(v)
as_tensor <- function(vector) {
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
  cbind(c(xx, yx, zx), c(yx, yy, zy), c(zx, zy, zz))
}
