#' Default regularization function for SSDiD
#' @description
#' The regularization term `\eta`^2 for SSDiD estimator proposed for usage in the SSDiD paper.
#'
#' You can create your own regularization function that uses s2 and N as arguments (s2 being the first parameter) and
#' pass it to sequential_estimator function.
#'
#' @param s2    Numeric. The noise variance.
#' @param N     Numeric.
#' @param deg   Numeric. The default value of 0.9 is proposed in the paper in Remark 3.2.
#'
#' @return      Numeric. The penalty term `eta^2` (regularization parameter).
#' @export
default_penalty <- function(s2, N) { s2 * N^(-0.9) }
