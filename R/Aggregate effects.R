#' Weight the estimated effect by cohort weight (population).
#' @param tau      Numeric matrix. The matrix with estimated effects.
#' @param W        Binary or boolen matrix. The matrix N x T of treatment indicators.
#' @param coh      Numeric vector. The N x 1 cohort weight (population) vector.
#' @param N0       Integer. The number of control (never treated) units.
#' @param N        Integer. The total number of units.
#' @export
aggregate_by_pop <- function(tau, W, coh, N0, N) {
  tau_avg <- array(0, dim = max(rowSums(W)))
  counter <- array(0, dim = max(rowSums(W)))

  for (t in 1:ncol(W)) {
    for (j in (N0 + 1):N) {
      if (W[j, t] == 1) {
        lag <- sum(W[j, 1:t])
        inc <- coh[j]
        tau_avg[lag] <- tau_avg[lag] + tau[j, t] * inc
        counter[lag] <- counter[lag] + inc
      }
    }
  }
  tau_lag <- tau_avg / counter
  tau_lag
}

#' Weight the estimated effect inversely proportional to DiD variance.
#' DiD variance is used instead of SSDiD variance to eliminate potential overfitting.
#' @param tau      Numeric matrix. The matrix with estimated effects.
#' @param W        Binary or boolen matrix. The matrix N x T of treatment indicators.
#' @param coh      Numeric vector. The N x 1 cohort weight (population) vector.
#' @param N0       Integer. The number of control (never treated) units.
#' @param N        Integer. The total number of units.
#' @export
# TODO: Should I use pi = share of popwt?
aggregate_inv_did_var <- function(tau, W, coh, N0, N) {
  tau_avg <- array(0, dim = max(rowSums(W)))
  counter <- array(0, dim = max(rowSums(W)))

  for (t in 1:ncol(W)) {
    for (j in (N0 + 1):N) {
      if (W[j, t] == 1) {
        lag <- sum(W[j, 1:t])
        inc <- 1 / ((1/sum(coh[1:(j-1)]) + 1/coh[j]) * (1 + 1/(t-1)))
        tau_avg[lag] <- tau_avg[lag] + tau[j, t] * inc
        counter[lag] <- counter[lag] + inc
      }
    }
  }
  tau_lag <- tau_avg / counter
  tau_lag
}
