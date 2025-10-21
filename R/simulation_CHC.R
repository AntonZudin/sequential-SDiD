#' Estimate the DGP for CHC simulation.
#' @description
#' The function estimates the DGP (interactive FE and noise) and adjusts the data to fit our needs:
#'   - Increase the number of units in each cohort
#'   - Adjust the ratio of noise variance to interactive FE variance
#' @importFrom dplyr select
#' @param panel        List.    The panel on the `unit` level.
#' @param signal_share Numeric. The share of interactive FE
#' @param scale        Integer. The number of units in each cohort is rescaled by this factor
#'                              to ensure our cohorts are large enough for the asymptotics to be relevant.
#' @param rank         Integer. The rank of interactive FE to extract from the data.
#' @return             List.    The list contains:
#'    - `L`        Dataframe. The dataframe contains interactive fixed effects.
#'    - `X`        Dataframe. The dataframe has auxiliary data like adoption date, population weights and covariates.
#'    - `sigma_2`  Numeric.   The noise variance on the unit level. The noise is supposed to be homoskedastic and does not have any autocorrelation.
#'    - `N`        Integer.   The total number of units after rescaling.
#'    - `T`        Integer.   The number of time periods.
#' @export estimate_dgp_chc
estimate_dgp_chc <- function(panel, signal_share = 0.0, scale = 4, rank = 5) {
  Y <- panel$Y; W <- panel$W; X <- panel$X
  N <- nrow(Y); T <- ncol(select(Y, -c("adopt_date")))

  model_m <- MCPanel::mcnnm_cv(M = as.matrix(select(Y, -c("adopt_date"))),
                    mask = 1 - as.matrix(select(W, -c("adopt_date"))))
  FE <- model_m$u %*% t(rep(1, T)) + rep(1, N) %*% t(model_m$v)
  L <- model_m$L
  {L_clean <- svd_compact(L)$u[, 1:rank] %*%
    svd_compact(L)$d[1:rank, 1:rank] %*% t(svd_compact(L)$v[, 1:rank])}
  eps_hat <- as.matrix(select(Y, -c("adopt_date")) - L_clean - FE)

  L_nrm <- L_clean/norm(L_clean, "F") * sqrt(N * T) * signal_share
  eps_nrm <- eps_hat/norm(eps_hat, "F") * sqrt(N * T) * (1 - signal_share)
  L_df <- data.frame(adopt_date = Y$adopt_date, L_nrm)
  #L_df <- as.data.frame(L_nrm)
  #L_df <- cbind(Y$adopt_date, L_df)
  colnames(L_df) <- colnames(Y); rownames(L_df) <- rownames(Y)

  sample_idx <- c()
  for (val in unique(Y$adopt_date)) {
    sample_idx <- c(sample_idx,
      sample(rownames(L_df[L_df$adopt_date == val,]),
             scale * nrow(L_df[L_df$adopt_date == val,]),
             replace = TRUE)
    )
  }
  L_sim <- L_df[sample_idx, ]
  X_sim <- X[sample_idx, ]
  W_sim <- W[sample_idx, ]
  sigma_2 <- var(as.vector(eps_nrm))
  list(L = L_sim, X = X_sim, W = W_sim,
       sigma_2 = sigma_2, N = nrow(L_sim), T = T)
}

#' Simulate the CHC DGP by sampling the noise.
#' @importFrom mvtnorm rmvnorm
#' @param params List.  The list should contain:
#'    - `L`        Dataframe. The dataframe contains interactive fixed effects.
#'    - `X`        Dataframe. The dataframe has auxiliary data like adoption date, population weights and covariates.
#'    - `sigma_2`  Numeric.   The noise variance on the unit level.
#'    - `N`        Integer.   The total number of units after rescaling.
#'    - `T`        Integer.   The number of time periods.
#' @return         List.      The list contains `Y` dataframe.
#' @export
simulate_dgp_chc <- function(params) {
  N <- params$N; T <- params$T
  sigma_2 <- params$sigma_2
  Y_sim <- params$L
  eps_sim <- mvtnorm::rmvnorm(n = N, mean = rep(0, T),
                              sigma = diag(sigma_2, nrow = T))
  # TODO Why does numbers for rownames work?
  {Y_sim[,names(Y_sim) %in% c(1959:1988)] <-
    Y_sim[,names(Y_sim) %in% c(1959:1988)] + eps_sim}
  list(Y = Y_sim)
}
