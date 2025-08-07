# FIXME: There is an issue that the aggregate_Y sorts the data by `adopt_date` in the increasing order,
#       so the correctness of the output depends on the right format of `adopt_date` column.
# TODO: Add sorting by colSums(W), so the unit order is correct

#' Aggregate dataframe on adoption_date and covariates.
#' @description
#' `NB`:  Dataframe Y should presorted in non increasing order by `adopt_date`.
#' @param Y         Dataframe. The dataframe should contain 'adopt_date', 'popwt', covariate and time columns with outcomes.
#' @param covs      Character vector. The covariate columns to be used in aggregation.
#' @param time_cols Character vector. The time columns of Y dataframe with outcomes.
#' @return          List. The list contains:
#'   - Y:  Dataframe. The dataframe has 'Group', 'popwt' and time columns. 'Group' column is concatenation of 'adopt_date' and covariate columns.
#'   - N:  Numeric. The number of rows of the Y dataframe after aggregation.
aggregate_Y <- function(Y, covs, time_cols) {

  group_vars <- c("adopt_date", covs)
  sum_cols <- c(time_cols, "popwt")

  if (length(group_vars) == 1) {
    group_factor <- Y[[group_vars[1]]]
  } else {
    group_factor <- interaction(Y[group_vars], sep = "_", drop = TRUE)
  }

  unique_groups <- levels(as.factor(group_factor))

  sum_data <- as.matrix(Y[sum_cols])
  aggregated_sums <- rowsum(sum_data, group_factor, na.rm = TRUE)

  group_info <- Y[match(unique_groups, group_factor), group_vars, drop = FALSE]

  Y_avg <- cbind(group_info, aggregated_sums)
  Y_avg$Group <- unique_groups

  if (length(covs) > 0) {
    Y_avg <- Y_avg[, !colnames(Y_avg) %in% covs, drop = FALSE]
  }
  N <- nrow(Y_avg)
  list(Y = Y_avg, N = N)
}

#' SVD decomposition of a matrix.
#' @param A Numeric matrix. SVD is applied on this matrix.
#' @return  List. The list contains:
#'   - `u`:  Numeric matrix. The left singular vectors of matrix A.
#'   - `d`:  Numeric vector. The singular values vector.
#'   - `v`:  Numeric matrix. The right singular vectors of matrix A
svd_compact <- function(A) {
  init <- svd(A)
  u <- init$u
  d <- init$d
  v <- init$v

  u_fin <- u[1:dim(A)[1], 1:qr(A)$rank]
  d_fin <- diag(d[1:qr(A)$rank], nrow = qr(A)$rank, ncol = qr(A)$rank)
  v_fin <- v[1:dim(A)[2], 1:qr(A)$rank]

  res <- list("u" = u_fin, "d" = d_fin, "v" = v_fin)
  res
}

# TODO: Write that we specify covariates by passing the arguments
# TODO: Write that if you want to avoid utilizing to_matrices function,
#       then you need to have "contr_cov_" and "treat_cov_"
# TODO: ?Sort unsorted dataframe, so that later treated units are later

#' Converts a data set into wide panels.
#' @param panel       Dataframe. The dataframe is utilized to
#' @param unit        Numeric or character. The column number or index that corresponds to the unit identifier.
#' @param time        Numeric or character. The column number or index that corresponds to the time identifier
#' @param outcome     Numeric or character. The column number or index that corresponds to the outcome identifier.
#' @param treatment   Numeric or character. The column number or index that corresponds to the treatment indicator.
#' @param population  Numeric, character or NULL. If not NULL, the column number or index that corresponds to the population identifier.
#' @param contr_covs  Numeric or character vector. If the length of the vector is not 0, the vector coresponds to the covariates for control (never-treated) units.
#' @param treat_covs  Numeric or character vector. If the length of the vector is not 0, the vector coresponds to the covariates for ever treated units.
#' @param never_treat Character. The time index that indicates never-treated units.
#' @return `panel`    List. The list contains:
#'   - `Y`:  Dataframe. This is a wide panel dataframe of outcomes with adoption date column being the first one.
#'   - `W`:  Dataframe. This is a wide panel dataframe of treatment indicators with adoption date column being the first one.
#'   - `X`:  Dataframe. This is a dataframe with auxiliary data like adoption date, population weights and covariates.
#' @export
to_matrices <- function(
    panel,
    unit = 1, time = 2, outcome = 3, treatment = 4, population = NULL,
    contr_covs = c(), treat_covs = c(), never_treat = "2500"
) {
  if (is.null(population)) {
    panel$popwt <- rep(1, nrow(panel))
    population  <- "popwt"
  }

  keep <- c(unit, time, outcome, treatment,
            population, contr_covs, treat_covs)

  contr_covs <- c(contr_covs); treat_covs <- c(treat_covs)

  if (!all(keep %in% 1:ncol(panel) | keep %in% colnames(panel))) {
    stop("Column identifiers should be either integer or column names in `panel`.")
  }

  index.to.name <- function(x) {
    new_names <- c()
    for (col in c(x)) {
      if (col %in% 1:ncol(panel)) {
        new_names <- c(new_names, colnames(panel)[col])
      } else {
        new_names <- c(new_names, col)
      }
    }
    if (length(new_names) == 1) {
      return(new_names[1])
    } else {
      return(new_names)
    }
  }

  unit <- index.to.name(unit)
  time <- index.to.name(time)
  outcome <- index.to.name(outcome)
  treatment <- index.to.name(treatment)
  population <- index.to.name(population)
  contr_covs <- unlist(sapply(contr_covs, index.to.name))
  treat_covs <- unlist(sapply(treat_covs, index.to.name))

  keep <- c(unit, time, outcome, treatment,
            population, contr_covs, treat_covs)

  panel <- panel[keep]
  if (!is.data.frame(panel)){
    stop("Unsupported input type `panel.`")
  }
  if (anyNA(panel)) {
    stop("Missing values in `panel`.")
  }
  if (length(unique(panel[, treatment])) == 1) {
    stop("There is no variation in treatment status.")
  }
  if (!all(panel[, treatment] %in% c(0, 1))) {
    stop("The treatment status should be in 0 or 1.")
  }
  # Convert potential factor/date columns to character
  panel <- data.frame(
    lapply(panel,
           function(col) {if (is.factor(col) || inherits(col, "Date")) as.character(col) else col}),
    stringsAsFactors = FALSE
  )
  val <- as.vector(table(panel[, unit], panel[, time]))
  if (!all(val == 1)) {
    stop("Input `panel` must be a balanced panel: it must have an observation for every unit at every time.")
  }

  panel <- panel[order(panel[, unit], panel[, time]), ]
  num.years <- length(unique(panel[, time]))
  num.units <- length(unique(panel[, unit]))

  ##unit level
  Y <- matrix(panel[, outcome], num.units, num.years, byrow = TRUE,
              dimnames = list(unique(panel[, unit]), unique(panel[, time])))
  W <- matrix(panel[, treatment], num.units, num.years, byrow = TRUE,
              dimnames = list(unique(panel[, unit]), unique(panel[, time])))
  unique_unit_id <- seq(1, num.units * num.years, num.years)

  X <- panel[unique_unit_id,
             c(population, contr_covs, treat_covs),
             drop = FALSE]
  names(X) <- c("popwt",
                vec_paste0("contr_cov_", length(contr_covs)),
                vec_paste0("treat_cov_", length(treat_covs)))
  row.names(X) <- unique(panel[, unit])

  unit_order <- order(rowSums(W))

  Y <- Y[unit_order, ]
  W <- W[unit_order, ]
  X <- X[unit_order, , drop = FALSE]

  #time_cols <- colnames(Y)
  treat_periods <- matrix(rowSums(W[, ]))
  periods <- colnames(W)

  adopt_date <- function(treat_period, periods, never_treat) {
    T <- length(periods)
    if (treat_period == 0) {
      return(never_treat)
    } else {
      return(periods[T + 1 - treat_period])
    }
  }
  ad_date <- apply(treat_periods, 1, adopt_date, periods, never_treat)

  Y <- data.frame(adopt_date = ad_date, Y, check.names = FALSE)
  W <- data.frame(adopt_date = ad_date, W, check.names = FALSE)
  X <- data.frame(adopt_date = ad_date, X, check.names = FALSE)

  #never_treat <- as.numeric(never_treat)

  list(Y =  Y, W = W, X = X,
       never_treat = never_treat)
}

# TODO: Check when Y and W can be a matrix
#' Aggregate panel to `cohort` level or prepare for the `unit` level.
#' @description
#' The function prepares the data for sequential_estimator.
#' For "unit" level the function removes adoption date column in Y and W. Y and W can be both matrices without `adopt_date` column.
#' For "cohort" level the function aggregates the outcome and population variable by adoption date and covariates.
#' @param panel List. The list should contain:
#' - `Y`    Dataframe or numeric matrix. This is a wide panel dataframe of outcomes with adoption date column being the first one. `Y` can be a matrix when the estimation is conducted on the `unit` level.
#' - `X`    Dataframe. This is a dataframe with auxiliary data like adoption date, population weights and covariates.
#' or
#' - `Y_wt` Dataframe. This is Y dataframe with merged X dataframe. It contains
#' - `W`    Dataframe, numeric matrix or NULL. This is a wide panel dataframe of treatment indicators with adoption date column being the first one. NULL can be passed when boot is TRUE. Numeric matrix can be passed when the estimation is conducted on the `unit` level.
#' @param level Character. The level should be `unit` or `cohort`. The data is not being aggregated when the level is "unit".
#' @param boot Bool. If TRUE, W_avg matrix is not calculated to speed up the bootstraping for standard error estimation.
#' @return `panel_avg` List. The list contains:
#'   -`Y_avg`:  Numeric matrix. The N x T matrix of outcomes.
#'   -`W_avg`:  Binary or boolen matrix. The N x T matrix of treatment indicators. The matrix has a stairlike structure with treated cells being in the bottom.
#'   -`coh`:    Numeric vector. The N x 1 cohort weights vector (the number of units in cohorts or the cohort population).
#' @export
prepare_matrices <- function(
    panel,
    level = "cohort",
    boot = FALSE
) {
  var_cond <- (all(c("Y", "X") %in% names(panel)) |
                 "Y_wt" %in% names(panel) &&
                 "never_treat" %in% names(panel) &&
                 ("W" %in% names(panel) | boot))
  if (!var_cond) {
    stop("Panel should contain `Y`, `X` or `Y_wt` and `W` dataframes and `never_treat` variable.")
  }

  if ("Y_wt" %in% names(panel)) {
    if (!("adopt_date" %in% names(panel$Y_wt))) {
      stop("`adopt_date` column should be present in Y_wt dataframe")
    }
  } else {
    col_cond <- (
        boot ||
        (level == "unit") ||
        (colnames(panel$W)[1] == "adopt_date" &&
           colnames(panel$Y)[1] == "adopt_date")
    )
    if (!col_cond) {
      stop("`adopt_date` should be the first column in both Y and W dataframes")
    }
  }

  if (level == "unit") {
    if ("Y" %in% names(panel) && "adopt_date" %in% colnames(panel$Y)) {
      panel$Y_avg <- as.matrix(subset(panel$Y, select = -adopt_date))
      if (!boot) {
        panel$W_avg <- as.matrix(subset(panel$W, select = -adopt_date))
      }
    } else if ("Y" %in% names(panel)) {
      panel$Y_avg <- panel$Y; if (!boot) panel$W_avg <- panel$W
    }
    panel$Y <- NULL; panel$W <- NULL
    if ("Y_wt" %in% names(panel)) panel$coh <- panel$Y_wt$popwt
    else panel$coh <- panel$X$popwt
    return(panel)
  }

  if (level == "unit") {
    if ("Y" %in% names(panel)) {
      if ("adopt_date" %in% colnames(panel$Y)) {
        panel$Y_avg <- as.matrix(subset(panel$Y, select = -adopt_date))
        if (!boot) {
          panel$W_avg <- as.matrix(subset(panel$W, select = -adopt_date))
        }
      } else {
        panel$Y_avg <- panel$Y; if (!boot) panel$W_avg <- panel$W
      }
      panel$coh <- panel$X$popwt
      panel$Y <- NULL; panel$W <- NULL
    } else {
      # Y_wt is utilized for Y_avg and coh
      Y_wt <- panel$Y_wt
      exclude_cols <- c("adopt_date", "popwt") # No covariates
      time_cols <- colnames(Y_wt)[!(colnames(Y_wt) %in% exclude_cols)]
      panel$Y_avg <- Y_wt[, time_cols]
      panel$coh <- panel$Y_wt$popwt
      if (!boot) {
        panel$W_avg <- as.matrix(subset(panel$W, select = -adopt_date))
      }
    }

    return(panel)
  }

  if (!("Y_wt" %in% names(panel))) {
    Y <- panel$Y
    X <- panel$X
    contr_covs <- names(X)[startsWith(names(X), "contr_cov_")]
    treat_covs <- names(X)[startsWith(names(X), "treat_cov_")]
    never_treat <- panel$never_treat

    time_cols <- colnames(Y)[-1]
    T <- length(time_cols)

    Y_wt <- merge(Y, X, by = c("row.names", "adopt_date"))
    row.names(Y_wt) <- Y_wt$`Row.names`
    Y_wt <- subset(Y_wt, select = -Row.names)
  } else {
    Y_wt <- panel$Y_wt
    contr_covs <- names(Y_wt)[startsWith(names(Y_wt), "contr_cov_")]
    treat_covs <- names(Y_wt)[startsWith(names(Y_wt), "treat_cov_")]
    never_treat <- panel$never_treat

    exclude_cols <- c("adopt_date", contr_covs, treat_covs, "popwt")
    time_cols <- colnames(Y_wt)[!(colnames(Y_wt) %in% exclude_cols)]
    T <- length(time_cols)
  }

  Y_wt[, time_cols] <- Y_wt[, time_cols] * Y_wt$popwt # paste weights
  contr_agg <- aggregate_Y(Y_wt[Y_wt$adopt_date == never_treat, ],
                           contr_covs,
                           time_cols)
  Y_avg_c <- contr_agg$Y
  N0 <- contr_agg$N

  treat_agg <- aggregate_Y(Y_wt[Y_wt$adopt_date != never_treat, ],
                           treat_covs,
                           time_cols)
  Y_avg_tr <- treat_agg$Y
  N1 <- treat_agg$N

  Y_avg <- rbind(Y_avg_tr, Y_avg_c)
  Y_avg[, time_cols] <- Y_avg[, time_cols] / Y_avg$popwt
  rownames(Y_avg) <- Y_avg$Group
  # TODO: Understand why we reverse everything
  coh <- rev(Y_avg$popwt)
  ad_date <- rev(Y_avg$adopt_date)
  Y_avg <- subset(Y_avg, select = -c(Group, popwt, adopt_date))
  Y_avg <- apply(Y_avg, 2, rev)


  if (!boot) {
    W_avg <- matrix(0, nrow = dim(Y_avg)[1], ncol = dim(Y_avg)[2],
                    dimnames = list(rownames(Y_avg), colnames(Y_avg)))
    for (i in (N0 + 1):(N1 + N0)) {
      t <- which(time_cols == ad_date[i])
      W_avg[i, t:T] <- 1
    }
  } else {
    W_avg <- NULL
  }

  return_panel <- list(
    Y_avg = Y_avg, W_avg = W_avg,
    coh = coh
  )

  if (!("Y_wt" %in% names(panel))) return_panel$X <- X
  if (!is.null(panel$N0)) return_panel$N0 <- panel$N0

  return_panel
}
