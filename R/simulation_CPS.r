#' Estimates the DGP parameters used in the placebo studies
#' of the Sequential Synthetic DiD paper.
#'
#' This placebo studies is an adjusted verion of Synthetic DiD paper.
#'
#' @param Y : Dataframe. It is a dataframe of outcomes.
#' @param assignment_vector : Binary/boolen vector. It is a Nx1 vector of treatment assignments.
#' @param rank : Integer. The rank of the estimated signal component L.
#'
#' @return A list with the following elements:
#'  - `F`:        Numeric matrix. The FE of the
#'  - `M`:        Numeric matrix. The interactive FE of the
#'  - `Sigma`:    Numeric matrix. The noise covariance matrix.
#'  - `sigma_2`: Numeric. The noise variance: the number which is located on the main diagonal of the covariance matrix.
#'  - `pi`:       Numeric vector.
#'  - `ar_coef`:  Numeric vector. AR(2) model coefficients for the covariance matrix Sigma.
#'        and an element ar_coef with the AR(2) model coefficients underlying the covariance Sigma
#'  - `rownames`: Character vector. The row names of `Y` matrix.
#'  - `colnames`: Character vector. The column names `Y` matrix.
#' @export estimate_dgp
# TODO: Do I need to save the column names in simulated Y as in the original Y?
estimate_dgp <- function(
  Y,
  assignment_vector,
  rank,
  pi_singular_vec,
  never_treat
) {
  N <- dim(Y)[1]
  T <- dim(Y)[2]
  # Convert Y to a matrix
  if (any(class(Y) == "data.frame")) {
    if ("adopt_date" %in% colnames(Y)) Y <- subset(Y, select = -c(adopt_date))
    Y <- as.matrix(Y)
  }
  overall_mean <- mean(Y)
  overall_sd <- norm(Y - overall_mean, 'f') / sqrt(N*T)
  Y_norm <- (Y - overall_mean)/overall_sd

  components <- decompose_Y(Y_norm, rank = rank)
  M <- components$M
  F <- components$F
  E <- components$E

  ar_coef <- round(fit_ar2(E), 2)
  cor_matrix <- ar2_correlation_matrix(ar_coef, T)
  scale_sd <- norm(t(E) %*% E /N, 'f') / norm(cor_matrix, 'f')
  cov_mat <- cor_matrix * scale_sd

  #treat_vector <- (assignment_vector != control_index)
  if (pi_singular_vec) {
    # In decompose_Y function unit_factors are multiplied by square root of the number of units.
    # We do min-max scaling, so the multiplication does not affect the result.
    pi <- decompose_Y(M, rank)$unit_factors[, 1]
  } else {
    pi <- glm(assignment_vector~unit_factors, family = 'binomial')$fitted.values
  }
  names(pi) <- rownames(Y)

  list(F = F, M = M, Sigma = cov_mat,
       sigma_2 = cov_mat[1, 1], never_treat = never_treat,
       pi = pi, ar_coef = ar_coef,
       colnames = colnames(Y), rownames = rownames(Y))
}


#' Simulates data from DGPs used in the placebo studies
#' of the synthetic difference in differences paper.
#' @importFrom mvtnorm rmvnorm
#'
#' @param parameters : List. The list contains dgp parameters (F, M, Sigma, pi) as output by estimate_dgp function.
#'                           F and M are numeric matrices.
#' @param N1 :         Integer. A cap on the number of treated units.
#' @param a_min :      Integer. The number of treated periods.
#' @param permute_pi : Bool. If TRUE, pi is permuted, so the treatment status is random.
#' @param as_df :      Bool. If TRUE, returns Y and W as dataframes instead of matrices.
#'
#' @return List with the following elements:
#'  - `Y`:       Numeric matrix or dataframe.
#'  - `W`:       Binary/boolen matrix or dataframe.
#'  - `Sigma`:   Numeric matrix.
#'  - `cohorts`: Numeric vector. The number of control units N0, and the number of control periods T0.
#'         The first N0 rows of Y are for units assigned to control, the remaining rows are for units assigned to treatment.
#' @export simulate_dgp
simulate_dgp <- function(params, N1, a_min, permute_pi = FALSE,  as_df = TRUE) {
  F <- params$F; M <- params$M
  Sigma <- params$Sigma
  pi <- params$pi
  if (permute_pi){
    pi <- sample(pi)
  }
  never_treat <- params$never_treat

  N <- nrow(M); T <- ncol(M)

  assignment <- randomize_treatment(pi, N, N1)

  if (permute_pi){
    pi <- sample(pi)
  }
  random_cohort <- randomize_cohort(
    assignment, pi,
    T, a_min
  )
  W <- random_cohort$W
  cohorts <- random_cohort$cohort

  N1 <- sum(assignment)

  Y <- F + M + mvtnorm::rmvnorm(N, sigma = Sigma)

  # Reorder the units in accordance with the treatment time:
  # The treatment indicator matrix (`W`) should have a stair-like structure
  # with treated cells being at the bottom.
  treatment_order <- order(cohorts, decreasing = TRUE)
  Y <- Y[treatment_order, ]
  W <- W[treatment_order, ]
  cohorts <- cohorts[treatment_order]

  rownames(Y) <- rownames(W) <- params$rownames
  colnames(Y) <- colnames(W) <- params$colnames
  if (as_df) {
    adopt_dates <- cohort_to_adopt_date(
      cohorts, params$colnames,
      never_treat, T
    )
    Y <- data.frame(adopt_date = adopt_dates, Y, check.names = FALSE)
    W <- data.frame(adopt_date = adopt_dates, W, check.names = FALSE)
  }
  list(
    Y = Y,
    W = W,
    Sigma = Sigma,
    cohort = cohorts
  )
}

#' Converts adoption dates from to column names of Y matrix.
#' @param cohort :       Numeric vector. The vector contains adoption dates
#' @param colnames_Y :   Character vector. The column names of Y matrix.
#' @param never_treat :  Character or Integer. The index that indicates control (never-treated) units.
#' @param T :            Integer. The number of time periods.
#'
#' @return `adopt_date`: Character or numeric vector.
cohort_to_adopt_date <- function(cohort, column_names, never_treat, T) {
  adopt_date <- cohort
  # Never-treated
  adopt_date[cohort > T] <- never_treat
  # Assign the adopt_date name according to the column names
  adopt_date[cohort <= T] <- column_names[cohort[cohort <= T]]
  adopt_date
}


#' Randomize treatment to n units with probability pi
#' then if the number of treated units is zero, assign treatment to one unit uniformly at random
#' and  if the number of treated units exceeds a cap, remove treatment uniformly at random so it is exactly that cap
#' @param pi : Numeric vector. The randomization probabilities.
#' @param N :  Integer. The total number of units.
#' @param N1 : Integer. The cap on the number of treated units.
#'
#' @return `assignment_sim` Binary vector. The vector is of length N, with ones indicating assignment to treatment.
randomize_treatment <-  function(pi, N, N1) {
  assignment_sim <- rbinom(N, 1, pi)
  index_as <- which(assignment_sim == 1)
  if (sum(assignment_sim) > N1) {
    index_pert <- sample(index_as, N1)
    assignment_sim <- rep(0, N)
    assignment_sim[index_pert] <- 1
  } else if (sum(assignment_sim) == 0){
    index_pert <- sample(1:N, N1)
    assignment_sim <- rep(0, N)
    assignment_sim[index_pert] <- 1
  }
  return(assignment_sim)
}


#' Randomize adoptions date for treated units.
#' @param assignment_vector : Binary vector. This is a binary vector of length N, with ones indicating assignment to treatment.
#' @param pi :                Numeric vector. The pi vector for units.
#' @param T :                 Integer. The total number of periods.
#' @param a_min :             Integer. The earliest possible treated cohort.
#'
#' @return List with the following elements:
#'  - `W`:      Binary or boolean matrix. It is a N x T treatment indicator matrix.
#'  - `cohort`: Integer vector. The assigned cohort vector.
randomize_cohort <- function(assignment_vector, pi, T, a_min) {

  N <- length(assignment_vector)

  # `T + 10` is never-treated cohort
  cohort <- rep(T + 10, N)
  scaled_pi <- (pi - min(pi)) / (max(pi) - min(pi))
  mu_i <- a_min + (T - a_min) * (1 - scaled_pi)

  #mu_i <- sample(mu_i) #reshuffle mu_i

  for (unit in 1:N) {
    if (assignment_vector[unit] != 0) {
      unit_cohort <- round(rnorm(n = 1, mean = mu_i[unit], sd = 1))
      unit_cohort <- max(unit_cohort, a_min) #It is not less than a_min
      unit_cohort <- min(unit_cohort, T) #Does not exceed T
      cohort[unit] <- unit_cohort
    }
  }

  W <- matrix(0, N, T)
  for (unit in 1:N) {
    if (cohort[unit] != T + 10) {
      W[unit, cohort[unit]:T] <- matrix(1, 1, (T - cohort[unit] + 1))
    }
  }

  return(list(W = W, cohort = cohort))
}


#' Decompose Y into components F, M, and E.
#' also computes a set of 'unit factors': the first [rank] left singular vectors from SVD(Y)
#' @param Y, the outcomes
#' @param rank, the assumed rank of the signal component L = F + M
#' @return a list with elements F, M, E, and unit_factors
decompose_Y <- function(Y, rank) {
  N <- dim(Y)[1]
  T <- dim(Y)[2]

  svd_data_mat <- svd(Y)
  factor_unit <- as.matrix(svd_data_mat$u[,1:rank]*sqrt(N))
  factor_time <- as.matrix(svd_data_mat$v[,1:rank]*sqrt(T))

  magnitude <- svd_data_mat$d[1:rank] / sqrt(N * T)
  L <- factor_unit %*% diag(magnitude, rank) %*% t(factor_time)

  E <- Y - L
  F <- outer(rowMeans(L), rep(1, T)) + outer(rep(1, N), colMeans(L)) - mean(L)
  M <- L - F

  return(list(F = F, M = M, E = E, unit_factors = factor_unit))
}

#' Estimate ar2 coefficients from iid time series
#' @param  E : Numeric matrix. It is a matrix of residuals with those time series as rows.

#' @return `ar_coef`: Numeric vector. It is a vector of ar2 coefficients: c(lag-1-coefficient, lag-2-coefficient).
fit_ar2 <- function(E){

  T_full <- dim(E)[2]
  E_ts <- E[, 3:T_full]
  E_lag_1 <- E[, 2:(T_full - 1)]
  E_lag_2 <- E[, 1:(T_full - 2)]

  a_1 <- sum(diag(E_lag_1 %*% t(E_lag_1)))
  a_2 <- sum(diag(E_lag_2 %*% t(E_lag_2)))
  a_3 <- sum(diag(E_lag_1 %*% t(E_lag_2)))

  matrix_factor <- rbind(c(a_1, a_3),c(a_3, a_2))

  b_1 <- sum(diag(E_lag_1 %*% t(E_ts)))
  b_2 <- sum(diag(E_lag_2 %*% t(E_ts)))

  ar_coef <- solve(matrix_factor) %*% c(b_1,b_2)
  return(ar_coef)
}

#' Compute the correlation matrix of a time series generated by an ar2 model.
#' @param ar_coef : Numeric vector. The coefficients of the ar2 model: c(lag-1-coefficient, lag-2-coefficient).
#' @param T :       Integer. The length of the time series.
#'
#' @return `cor_matrix` : Numeric matrix. The correlation matrix.
ar2_correlation_matrix <- function(ar_coef, T) {

  result <- rep(0, T)
  result[1] <- 1
  result[2] <- ar_coef[1]/(1 - ar_coef[2])
  for (t in 3:T){
    result[t] <- ar_coef[1]*result[t - 1] + ar_coef[2]*result[t - 2]
  }

  index_matrix <- outer(1:T, 1:T, function(x, y){ abs(y - x) + 1})
  cor_matrix <- matrix(result[index_matrix], ncol = T, nrow = T)

  return(cor_matrix)
}
