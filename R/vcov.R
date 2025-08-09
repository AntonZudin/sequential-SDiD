#' Compute standard errors for the sequential estimator.
#' @description
#' The standard errors are computed with Bayesian bootstrap. The weights are utilized in population weights.
#' @param object       `sequential_estimate` class object.
#' @param panel        List. The list should contain:
#'   - `Y`:  Dataframe. This is a wide panel dataframe of outcomes with adoption date column being the first one.
#'   - `W`:  Dataframe. This is a wide panel dataframe of treatment indicators with adoption date column being the first one.
#'   - `X`:  Dataframe. This is a dataframe with auxiliary data like adoption date, population weights and covariates.
#'   - `never_treat`:  Character. The index of control (never-treated) cohort.
#' @param B            Integer. The number of bootstrap replications.
#' @param type         Character or NULL. The type should be `did`, `sdid` or `both`.
#' @param return_tau_b Bool. If TRUE, returns the bootstrap replications.
#' @return             Numeric vector or List.
#' @export
# TODO: Understand how to use ... in R to pass extra args to function.
# TODO: Think what should I do with the name vcov while returning standard error
vcov.sequential_estimate <- function(
  object,
  panel,
  B = 1000,
  type = NULL,
  return_tau_b = FALSE
) {
  if (is.null(type)) { type <- attr(object, "type") }
  level <- attr(object, "level")
  s2 <- attr(object, "s2")

  if (!level %in% c("cohort", "unit")) {
    stop("Level should be either `cohort` or `unit`.")
  }
  if (is.null(s2)) {
    stop("s2 is not present as an attribute of `sequential_estimate` class .")
  }

  if (!(type %in% c("did", "sdid", "both"))){
    stop("Type should be either `did`, `sdid` or `both`.")
  }

  Y <- panel$Y
  # Remove adoption date from X dataframe
  if ("adopt_date" %in% colnames(panel$X)) {
    X <- subset(panel$X, select = -adopt_date)
  } else { X <- panel$X}

  both_est <- (type == "both")
  N <- dim(panel$Y)[1]
  W_avg <- prepare_wide(panel, level = level)$W_avg

  # TODO: Rename variables to make them more intuitive
  tau_lag_b <- array(dim = c(max(rowSums(W_avg[, -1])), B))
  if (both_est) tau_lag_b_2 <- array(dim = c(max(rowSums(W_avg[, -1])), B))

  N0 <- sum(rowSums(W_avg) == 0)
  panel$N0 <- N0
  panel$W <- NULL
  Y_wt <- merge(Y, X, by = "row.names", sort = FALSE)

  rownames(Y_wt) <- Y_wt$`Row.names`
  panel$Y_wt <- subset(Y_wt, select = -Row.names)
  panel$Y <- panel$X <- NULL

  for (b in 1:B) {
    wts_g <- rexp(n = N, rate = 1)
    wts <- wts_g/sum(wts_g) * N
    panel_b <- panel

    panel_b$Y_wt$popwt <- panel_b$Y_wt$popwt * wts

    panel_b_agg <- prepare_wide(panel_b, level, TRUE)
    panel_b_agg$W_avg <- W_avg
    # TODO: Should panel be replaced with more explicit Y, W and coh?
    estimate <- sequential_estimator(
      panel_b_agg, level, s2, type
    )
    if (both_est) {
      tau_lag_b[, b] <- estimate[["tau_sdid"]]
      tau_lag_b_2[, b] <- estimate[["tau_did"]]
    } else {
      tau_lag_b[, b] <- estimate
    }
  }
  result <- apply(tau_lag_b, 1, sd)
  if (both_est) {
    tau_se_2 <- apply(tau_lag_b_2, 1, sd)
    result <- list(se_sdid = result, se_did = tau_se_2)
  }

  if (return_tau_b && both_est) {
    result <- list(
      se = result,
      tau_b = list(tau_b_sdid = tau_lag_b, tau_b_did = tau_lag_b_2)
    )
  } else if (return_tau_b && !both_est) {
    result <- list(se = result, tau_b = tau_lag_b)
  }
  return(result)

  # TODO: Decide how to implement the return of se better
  # Add standard errors
  #se_est <- apply(tau_lag_b, 1, sd)
  #if (both_est) {
  #  se_did <- apply(tau_lag_b_2, 1, sd)
  #  se_est <- list(
  #    se_sdid = se_est,
  #    se_did = se_did
  #  )
  #}
  #attr(object, "se") <- se_est

  # Add bootstrapped tau
  #if (return_tau_b) {
  #  tau_b <- tau_lag_b
  #  if (both_est) {
  #    tau_b <- list(
  #      tau_b_sdid = tau_b,
  #      tau_b_did = tau_lag_b_2
  #    )
  #  }
  #  attr(object, "tau_b") <- tau_b
  #}

  #return(object)
}
