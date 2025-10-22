#' The list with estimator names.
#' It is created for convenience.
est_names <- list(
  sdid = "Synthetic DiD",
  did = "DiD"
)

#' Print a sequential_estimate object
#'
#' @param x : The sequential_estimate object to print.
#' @param digits : Integer. The number of digits to print after the dot.
#' @param n_lags : Integer. The number of lags to print out.
#' @param ... : Additional arguments (for compatibility with generic `print`).
#' @method print sequential_estimate
#' @export
print.sequential_estimate <- function(x, digits = 2, n_lags = NULL, ...) {
  type <- attr(x, "type")
  if (is.null(n_lags)){
    n_lags <- switch(type,
                both = length(x[[1]]),
                length(x))
  }
  fmt <- paste0("%1.", digits, "f")

  if (type == "both") {
    cat("Synthetic DiD:\n")
    cat(sprintf(fmt, x$tau_sdid[1:n_lags]), "\n")

    cat("\nDiD:\n")
    cat(sprintf(fmt, x$tau_did[1:n_lags]), "\n")
  } else {
    cat(est_names[[type]], ":\n", sep = "")
    cat(sprintf(fmt, unlist(x)[1:n_lags]), "\n")
  }

  # Print s2
  cat("\ns2:", sprintf(fmt, attr(x, "s2")), "\n")

  invisible(x)
}
