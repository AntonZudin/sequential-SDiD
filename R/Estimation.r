#' Penalty for SSDiD
#' @description
#' The regularization term `\eta`^2 for SSDiD estimator without s2.
#' @param N     Numeric.
#' @param deg   Numeric. The default value of 0.9 is proposed in the paper in Remark 3.2.
#' @return      Numeric. The penalty term without s2.
penalty <- function(N, deg = 0.9) N^(-deg)

#' Remove `adopt_date` column from a dataframe (Y or W).
#' @description
#' If the argument is a matrix, none of the columns are dropped.
#' @param df Dataframe/matrix. Drop `adopt_date` from this dataframe.
remove_adopt_date <- function(df) {
  if (is.null(colnames(df))) {
    return(df)
  } else {
    return(
      df[, !(colnames(df) %in% "adopt_date")]
    )
  }
}

#' Take the last column of a dataframe or matrix.
#' @param X  Dataframe or matrix.
#' @export
last_col = function(X) {
  if (!any(class(X) %in% c("data.frame", "matrix"))) {
    stop("X should be either data.frame or matrix.")
  }
  X[, ncol(X)]
}


#' Concatenates the `base_string` with the numbers from 1 to `len`.
#' @description
#' Creates a character vector containing strings of the following type "{base_string}{i}",
#' where `i` is a number between 1 and `len`.
#' If `len` argument is not positive, a vector of zero length is returned.
#' @param base_string : Character. The base character that is used in concatenation.
#' @param len :         Integer. The length of the output vector.
#'
#' @return              Character vector.
vec_paste0 <- function(base_string, len) {
  if (len > 0) {
    return(paste0(base_string, 1:len))
  } else {
    return(c())
  }
}


#' Computes the base synthetic diff-in-diff or diff-in-diff estimate.
#' @description
#' The bottom right cell is the only cell being treated (W_it = 1).
#' @param Y :    Numeric matrix. A submatrix of outcomes with one treated obs in the bottom right corner.
#' @param n_j :  Numeric vector. The vector of cohort weights (the number of units in cohorts or the cohort population).
#' @param s2 :   Numeric. The upper bound estimate of noise variance.
#' @param type : Character. Type of the estimator should be `sdid` or `did`.
#'
#' @return `tau` : Numeric. The base estimate on the Y submatrix.
#'
tau_sdid <- function(Y, n_j, s2, type = "sdid") {
  if (nrow(Y) != length(n_j)) {
    stop("The number of cohorts does not coincide in Y matrix and n_j vector")
  }
  if (!(type %in% c('did', 'sdid'))) {
    stop("The 'type' argument should be either 'sdid' or 'did'")
  }

  j_c <- nrow(Y) - 1
  t_c <- ncol(Y) - 1
  N <- sum(n_j[1:j_c])
  pi <- n_j[1:j_c] / N

  if (j_c > 1) {
    Y_c <- Y[1:j_c, 1:t_c]
  } else {
    Y_c <- t(Y[1:j_c, 1:t_c])
  }

  Y_j0 <- Y[j_c + 1, 1:t_c]
  Y_t0 <- Y[1:j_c, t_c + 1]
  Y_j0_t0 <- Y[j_c + 1, t_c + 1]

  if (type == "did") {
    lambda_reg <- rep(1/t_c, t_c)
    gamma_reg <- pi
  } else {
    # Precompute variables that are used multiple times
    pen <- penalty(N)
    ones_t <- rep(1, t_c)
    ones_j <- rep(1, j_c)

    Sigma_tc <- diag(s2 / pi, nrow = j_c)

    grad_gamma1 <- 2 * (Y_c %*% (-Y_j0))
    grad_gamma2 <- 2 * sum(-Y_j0)
    grad_reg <- c(grad_gamma1, grad_gamma2, -1)

    block_1 <- 2*tcrossprod(Y_c) + 2*pen * Sigma_tc
    block_2 <- 2*(Y_c %*% ones_t)
    hess_reg <- rbind(
      cbind(block_1, block_2, ones_j),
      c(block_2, 2 * t_c, 0),
      c(ones_j, 0, 0)
    )
    # This prevents the code from breaking,
    # forcing to do DID when hessian of the objective function is singular
    if (rcond(hess_reg) < 2.5e-16) gamma_reg <- pi
    else gamma_reg <- solve(hess_reg, -grad_reg)[1:j_c]

    diag_val_jc <- (1/j_c) * s2 * sum(1/pi)
    Sigma_jc <- diag(diag_val_jc, nrow = t_c)

    grad_lambda_1 <- 2*crossprod(Y_c, -Y_t0)
    grad_lambda_2 <- 2*sum(-Y_t0)
    gradl_reg <- c(grad_lambda_1, grad_lambda_2, -1)

    blockl_1 <- 2*crossprod(Y_c) + 2*pen * Sigma_jc
    blockl_2 <- 2*(crossprod(Y_c, ones_j))
    hessl_reg <- rbind(
      cbind(blockl_1, blockl_2, ones_t),
      c(blockl_2, 2*j_c, 0),
      c(ones_t, 0, 0)
    )

    if (rcond(hessl_reg) < 2.5e-16) lambda_reg <- rep(1/t_c, t_c)
    else lambda_reg <- solve(hessl_reg, -gradl_reg)[1:t_c]
  }

  tau <- (Y_j0_t0 - crossprod(Y_j0, lambda_reg)[1, 1]) -
    (crossprod(Y_t0, gamma_reg)[1, 1] - (t(gamma_reg) %*% Y_c %*% lambda_reg)[1, 1])
  tau
}

#' Computes the `cohort` level sequential synthetic diff-in-diff or diff-in-diff estimate.
#' @description
#' The outcome variable is aggregated on cohort or covariate level.
#' @param Y_avg :     Numeric matrix. The NxT aggregated matrix of outcomes.
#' @param W_avg :     Binary or boolean matrix. The N x T matrix of treatment indicators.
#' @param coh :       Numeric vector. The Nx1 cohort weights vector (the number of units in cohorts or the cohort population).
#' @param s2 :        Numeric. The upper bound estimate of noise variance.
#' @param type :      Character. Type of the estimator should be `sdid` or `did`.
#' @param N0 :        Integer. The number of control (never-treated) units.
#'
#' @return `tau_lag`: Numeric vector. The max_lag x 1 vector of treatment effects aggregated across units.
#'
estimation_cohort <- function(Y_avg, W_avg, coh, s2, type, N0 = 1) {
  if (!(type %in% c('did', 'sdid'))) {
	  stop("The 'type' argument should be either 'sdid' or 'did'")
  }

  tau_hat <- matrix(0, nrow = nrow(W_avg), ncol = ncol(W_avg))

  for (t in 1:ncol(W_avg)) {
    for (j in (N0 + 1):nrow(W_avg)) {
      if (W_avg[j, t] == 1) {
        tau_est <- base_estimator(Y_avg[1:j, 1:t], coh[1:j],
                            s2 = s2, type = type)
        tau_hat[j, t] <- tau_est
        Y_avg[j, t] <- Y_avg[j, t] - tau_est
      }
    }
  }
  tau_hat
}


#' Computes the `unit` level sequential synthetic diff-in-diff or diff-in-diff estimate.
#' @description
#' The outcome variable is not aggregated or the data was initially preaggregated on a high level.
#' @importFrom dplyr %>% group_by summarize arrange row_number desc
#' @param Y_avg :     Numeric matrix. The NxT matrix of outcomes.
#' @param W_avg :     Binary or boolean matrix. The NxT matrix of treatment indicators.
#' @param coh :       Numeric vector. The Nx1 cohort weights vector (the number of units in cohorts or the cohort population).
#' @param s2 :        Numeric. The upper bound estimate of noise variance.
#' @param type :      Character. Type of the estimator should be `sdid` or `did`.
#' @param N0 :        Integer. The number of control (never-treated) units.
#'
#' @return `tau_lag`: Numeric vector. The max_lag x 1 vector of treatment effects aggregated across units.
#'
estimation_unit <- function(Y, W, coh, s2, type, N0 = 1) {
  # TODO: Rename the variables to comprehend better their purpose
  if (!(type %in% c('did', 'sdid'))) {
	  stop("The 'type' argument should be either 'sdid' or 'did'")
  }

  N <- nrow(Y); T <- ncol(Y)
  tau_hat <- matrix(0, nrow = N, ncol = T)
  adopt_date <- T + 1 - rowSums(W)
  # TODO: Fix the bug that `dplyr` does not load
  cohort_df <- data.frame(adopt_date = adopt_date) %>%
    group_by(adopt_date) %>%
    summarize(n = dplyr::n()) %>%
    arrange(desc(row_number()))

  for (t in 1:T) {
    treat_mask <- cohort_df$adopt_date <= t
    n_coh_tr <- sum(treat_mask)
    if (n_coh_tr == 0) next

    cohorts <- cohort_df[treat_mask, ]
    j_c <- N - sum(cohorts$n)
    cont_idx <- 1:j_c
    # Process each cohort in descending order
    for (c in 1:n_coh_tr) {
      size_tr_coh <- cohorts$n[[c]]
      for (j in 1:size_tr_coh) {
        idx <- c(cont_idx, j_c + j)
        tau_est <- base_estimator(Y[idx, 1:t], coh[idx],
                            s2 = s2, type = type)
        tau_hat[j_c + j, t] <- tau_est
        Y[j_c + j, t] <- Y[j_c + j, t] - tau_est
      }
      j_c <- j_c + size_tr_coh
    }
  }
  tau_hat
}

#' `unit` and `cohort` level estimation functions.
#' @description
#' List containing `estimation_cohort` and `estimation_few` functions.
#' The first element is `estimation_cohort` and the second one is `estimation_few`.
#' This list is created for convenience.
estimation_funcs <- list(
  cohort = estimation_cohort,
  unit = estimation_unit
)


#' Computes the `unit` or `cohort` level sequential synthetic diff-in-diff or diff-in-diff estimate.
#' @description
#' Basically, this function is a wrapper which calls either `estimation_cohort` or `estimation_unit` function.
#' Then it aggregates the
#' @param panel_avg : List. The list should contain:
#'   - `Y_avg`: Numeric matrix. The N x T matrix of outcomes.
#'   - `W_avg`: Binary or boolen matrix. The N x T matrix of treatment indicators. The matrix has a stair-like structure with treated cells being at the bottom.
#'   - `coh`:   Numeric vector. The Nx1 cohort weights vector (the number of units in cohorts or the cohort population).
#'   - `X`:     Dataframe. The dataframe is not aggregated to cohort and covariate level. It contains auxiliary data like adoption date, population weights and covariates.
#' @param level :            Character. The level should be `unit` or `cohort`.
#' @param s2 :               Numeric or NULL. The upper bound estimate of noise variance.
#' @param type :             Character. Type of the estimator should be `sdid`, `did` or `both`.
#' @param return_s2 :        Bool. If TRUE, returns s2. If s2 is not passed (NULL), s2 is estimated and returned.
#' @param aggregate_effect : Function. The function aggregates the treatment effect for every lag by weighting the effect for every cohort.
#'                           The function should contain the following arguments:
#'   - `tau`: Numeric matrix or array.
#'   - `W`:   Numeric matrix.
#'   - `coh`: Numeric vector.
#'   - `N0`:  Integer. The number of control (never-treated) units.
#'   - `N`:   Integer. The total number of units.
#'
#' @return                 sequential_estimator
#' @export
sequential_estimator <- function(
  panel_avg,
  level,
  s2 = NULL,
  type = "sdid",
  aggregate_effect = aggregate_inv_did_var
) {
  if (!(level %in% c("cohort", "unit"))) {
    stop("The estimation level should be either 'cohort' or 'unit'.")
  }

  if (!(type %in% c("did", "sdid", "both"))) {
    stop("The estimation type should be either 'did','sdid' or 'both'.")
  }

  Y <- panel_avg$Y_avg; W <- panel_avg$W_avg
  N <- nrow(panel_avg$Y_avg)
  coh <- panel_avg$coh

  # TODO: Check whether N0 calculation should be removed
  if (!is.null(panel_avg$N0)) {
    N0 <- panel_avg$N0
  } else {
    N0 <- sum(rowSums(W) == 0)
  }

  if (is.null(s2)) s2 <- estimate_s2(Y, W)

  if (type == "both") {
    tau_hat_sdid <- estimation_funcs[[level]](Y, W, coh,
                      s2, "sdid", N0)
    tau_hat_did <- estimation_funcs[[level]](Y, W, coh,
                    s2, "did", N0)
    tau_sdid <- aggregate_effect(tau = tau_hat_sdid, W = W, coh = coh,
                                 N0 = N0, N = N)
    tau_did <- aggregate_effect(tau = tau_hat_did, W = W, coh = coh,
                                N0 = N0, N = N)
    estimate <- list(tau_sdid = tau_sdid, tau_did = tau_did)
  } else {
    # TODO: Change to correct version
    tau_hat <- estimation_funcs[[level]](Y, W, coh,
                  s2, type, N0)
    estimate <- aggregate_effect(tau = tau_hat, W = W, coh = coh,
                                 N0 = N0, N = N)
  }

  class(estimate) <- "sequential_estimate"
  attr(estimate, "type") <- type
  attr(estimate, "level") <- level
  attr(estimate, "s2") <- s2
  estimate
}


#' Estimates the upper bound of the noise variance.
#' @description
#' Removes time and unit fixed effects and computes the variance of residuals.
#'
#' There are 3 estimate options:
#' 1. OLS without population weigths. Do not pass `pop` vector and  set `wls` to FALSE.
#' 2. OLS with population weights. Pass `pop` vector and  set `wls` to FALSE.
#' 3. WLS utilizing population weigths. Pass `pop` vector and set `wls` to TRUE.
#'
#' @param Y Dataframe or matrix. The disaggregated wide panel of outcomes.
#' @param W Dataframe or matrix. The disaggregated wide panel of treatment indicators.
#' @param pop Numeric vector or NULL. The population weight vector.
#' @param wls Bool. If TRUE, WLS is used to to estimate TWFE.
#' @return `s2` Numeric. The estimated upper bound of the noise variance.
#' @export
#'
# TODO: Should I add an explicit argument for weighting
estimate_s2 <- function(Y, W, pop = NULL, wls = FALSE) {
  if (!requireNamespace("fixest", quietly = TRUE)) {
      stop("Please install the 'fixest' package: install.packages('fixest')")
  }
  Y <- remove_adopt_date(Y); W <- remove_adopt_date(W)

  wls <- !is.null(pop) && wls
  weight <- !is.null(pop) && !wls

  n_units <- nrow(Y)
  n_periods <- ncol(Y)

  if (any(class(Y) == "data.frame") && any(class(W) == "data.frame")) {
    Y <- unlist(Y); W <- unlist(W)
  } else if (!(any(class(Y) == "matrix") && any(class(W) == "matrix"))) {
    stop("`Y` and `W` should both be either matrices or dataframes.")
  }

  data_long <- data.frame(
    Y = as.vector(Y),
    unit = factor(rep(1:n_units, times = n_periods)),
    time = factor(rep(1:n_periods, each = n_units)),
    W = as.vector(W)
  )
  if (!is.null(pop)) data_long$pop <- rep(pop, times = n_periods)

  data_untreated <- subset(data_long, W == 0)
  if (nrow(data_untreated) == 0) {
    stop("No untreated observations found.")
  }

  n <- nrow(data_untreated)
  if (!is.null(pop)) {
    data_untreated$weights <- (data_untreated$pop / sum(data_untreated$pop)) * n
  }

  n_units <- nlevels(data_untreated$unit)
  n_periods <- nlevels(data_untreated$time)
  deg_free <- n - n_units - n_periods

  if (deg_free <= 0) {
    stop("Insufficient data: Degrees of freedom <= 0.")
  }

  # Remove fixed effects
  if (!wls) {
    model <- fixest::feols(Y ~ 0 | unit + time,
      data = data_untreated
    )
    if (weight) s2 <- sum(model$residuals^2 * data_untreated$weights) / deg_free
    else s2 <- model$ssr / deg_free

  } else {
    model <- fixest::feols(Y ~ 0 | unit + time,
      data = data_untreated,
      weights = ~weights
    )
    s2 <- model$ssr / deg_free
  }
  s2
}
