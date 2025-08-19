#' The list with estimator names.
#' It is created for convenience.
est_names <- list(
  sdid = "Synthetic DiD",
  did = "DiD"
)

#' Print a sequential_estimate object
#' @param x : The sequential_estimate object to print.
#' @param digits : Integer. The number of digits to print after the dot.
#' @param ... : Additional arguments (for compatibility with generic `print`).
#' @method print sequential_estimate
#' @export
print.sequential_estimate <- function(x, digits = 2, ...) {
  type <- attr(x, "type")
  fmt <- paste0("%1.", digits, "f")

  if (type == "both") {
    cat("Synthetic DiD:\n")
    cat(sprintf(fmt, x$tau_sdid), "\n")

    cat("\nDiD:\n")
    cat(sprintf(fmt, x$tau_did), "\n")
  } else {
    cat(est_names[[type]], ":\n", sep = "")
    cat(sprintf(fmt, unlist(x)), "\n")
  }

  # Print s2
  cat("\ns2:", sprintf(fmt, attr(x, "s2")), "\n")

  invisible(x)
}
