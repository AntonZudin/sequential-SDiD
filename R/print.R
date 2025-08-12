#' The list with estimator names.
#' It is created for convenience.
est_names <- list(
  sdid = "Synthetic DiD",
  did = "DiD"
)


#' Print a sequential_estimate object
#' @param x : The sequential_estimate object to print.
#' @param ... Additional arguments (for compatibility with generic `print`).
#' @method print sequential_estimate
#' @export
print.sequential_estimate <- function(x, ...) {
  type <- attr(x, "type")

  if (type == "both") {
    cat("Synthetic DiD:\n")
    print(x$tau_sdid)

    cat("\nDiD:\n")
    print(x$tau_did)
  } else {
    cat(est_names[[type]], ":\n", sep = "")
    print(x[1:length(x)])
  }

  cat("\ns2:", attr(x, "s2"), "\n")

  invisible(x)
}
