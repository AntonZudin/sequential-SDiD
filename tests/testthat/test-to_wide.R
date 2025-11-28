permute_unit <- function(unit_ids, T) {
  permut <- rep(0, length(unit_ids) * T)
  for (i in 1:length(unit_ids)){
    permut[(T*(i - 1) + 1):(T*i)] <- (T*(unit_ids[i] - 1) + 1):(T*unit_ids[i])
  }
  permut
}

permute_time <- function(time_ids, N) {
  permut <- rep(0, length(time_ids) * N)
  T <- length(time_ids)
  for (i in 1:N){
    permut[(T*(i - 1) + 1):(T*i)] <- T*(i - 1) + time_ids
  }
  permut
}


test_that("From long to wide works fine without population, without covariates", {
  panel_long <- data.frame(
    outcome = 1:100,
    unit = rep(1:10, each = 10),
    time = rep(1:10, 10),
    w = c(rep(0, 60),
          rep(0, 6), rep(1, 4),
          rep(c(rep(0, 4), rep(1, 6)), 2),
          rep(0, 2), rep(1, 8))
  )

  # Permute the ever-treated rows
  permut <- c(1, 2, 3, 4, 5, 6, 9, 7, 10, 8)
  panel_long <- panel_long[permute_unit(permut, 10), ]

  # Permute the columns
  permut <- c(10, 2, 3, 4, 5, 6, 7, 8, 9, 1)
  panel_long <- panel_long[permute_time(permut, 10), ]

  panel <- to_wide(panel_long, unit = "unit", time = "time", outcome = 1, treatment = 4, sort = TRUE)

  adopt_date <- c(rep("2500", 6), "7", "5", "5", "3")

  Y <- matrix(1:100, byrow = TRUE, nrow = 10)
  Y_correct <- data.frame(
    adopt_date = adopt_date,
    Y,
    check.names = FALSE
  )
  rownames(Y_correct) <- as.character(1:10)

  W <- matrix(
    c(rep(0, 60),
      rep(0, 6), rep(1, 4),
      rep(c(rep(0, 4), rep(1, 6)), 2),
      rep(0, 2), rep(1, 8)),
    byrow = TRUE,
    nrow = 10
  )
  W_correct <- data.frame(
    adopt_date = adopt_date,
    W,
    check.names = FALSE
  )
  rownames(W_correct) <- as.character(1:10)
  X_correct <- data.frame(adopt_date = adopt_date, popwt = rep(1, 10))

  expect_s3_class(panel$Y, "data.frame")
  expect_s3_class(panel$W, "data.frame")
  expect_s3_class(panel$X, "data.frame")
  expect_equal(panel$Y, Y_correct)
  expect_equal(panel$W, W_correct)
  expect_equal(panel$X, X_correct)
})


test_that("From long to wide works fine with population, without covariates", {
  panel_long <- data.frame(
    outcome = 1:100,
    unit = rep(1:10, each = 10),
    time = rep(1:10, 10),
    w = c(rep(0, 60),
          rep(0, 6), rep(1, 4),
          rep(c(rep(0, 4), rep(1, 6)), 2),
          rep(0, 2), rep(1, 8)),
    population = rep(21:30, each = 10)
  )

  # Permute the ever-treated rows
  permut <- c(1, 2, 3, 4, 5, 6, 9, 7, 10, 8)
  panel_long <- panel_long[permute_unit(permut, 10), ]

  panel_long[order(panel_long[, "unit"], panel_long[, "time"]), ]

  # Permute the columns
  permut <- c(10, 2, 3, 4, 5, 6, 7, 8, 9, 1)
  panel_long <- panel_long[permute_time(permut, 10), ]

  panel <- to_wide(
    panel_long,
    unit = "unit", time = "time",
    outcome = 1, treatment = 4,
    weights = 5,
    sort = TRUE)

  adopt_date <- c(rep("2500", 6), "7", "5", "5", "3")

  Y <- matrix(1:100, byrow = TRUE, nrow = 10)
  Y_correct <- data.frame(
    adopt_date = adopt_date,
    Y,
    check.names = FALSE
  )
  rownames(Y_correct) <- as.character(1:10)

  W <- matrix(
    c(rep(0, 60),
      rep(0, 6), rep(1, 4),
      rep(c(rep(0, 4), rep(1, 6)), 2),
      rep(0, 2), rep(1, 8)),
    byrow = TRUE,
    nrow = 10
  )
  W_correct <- data.frame(
    adopt_date = adopt_date,
    W,
    check.names = FALSE
  )
  rownames(W_correct) <- as.character(1:10)
  X_correct <- data.frame(adopt_date = adopt_date, popwt = 21:30)

  expect_s3_class(panel$Y, "data.frame")
  expect_s3_class(panel$W, "data.frame")
  expect_s3_class(panel$X, "data.frame")
  expect_equal(panel$Y, Y_correct)
  expect_equal(panel$W, W_correct)
  expect_equal(panel$X, X_correct)
})


test_that("From long to wide works fine with population, with covariates", {
  panel_long <- data.frame(
    outcome = 1:100,
    unit = rep(1:10, each = 10),
    time = rep(1:10, 10),
    w = c(rep(0, 60),
          rep(0, 6), rep(1, 4),
          rep(c(rep(0, 4), rep(1, 6)), 2),
          rep(0, 2), rep(1, 8)),
    population = rep(21:30, each = 10),
    contr_cov_1 = c(rep(1, 30), rep(2, 30), rep(1, 40)),
    treat_cov_1 = c(rep(1, 60), rep(1, 10), rep(2, 20), rep(1, 10))
  )

  # Permute the ever-treated rows
  permut <- c(1, 2, 3, 4, 5, 6, 9, 7, 10, 8)
  panel_long <- panel_long[permute_unit(permut, 10), ]

  panel_long[order(panel_long[, "unit"], panel_long[, "time"]), ]

  # Permute the columns
  permut <- c(10, 2, 3, 4, 5, 6, 7, 8, 9, 1)
  panel_long <- panel_long[permute_time(permut, 10), ]

  panel <- to_wide(
    panel_long,
    unit = "unit", time = "time",
    outcome = 1, treatment = 4,
    weights = "population",
    contr_covs = "contr_cov_1",
    treat_covs = c("treat_cov_1"),
    sort = TRUE)

  adopt_date <- c(rep("2500", 6), "7", "5", "5", "3")

  Y <- matrix(1:100, byrow = TRUE, nrow = 10)
  Y_correct <- data.frame(
    adopt_date = adopt_date,
    Y,
    check.names = FALSE
  )
  rownames(Y_correct) <- as.character(1:10)

  W <- matrix(
    c(rep(0, 60),
      rep(0, 6), rep(1, 4),
      rep(c(rep(0, 4), rep(1, 6)), 2),
      rep(0, 2), rep(1, 8)),
    byrow = TRUE,
    nrow = 10
  )
  W_correct <- data.frame(
    adopt_date = adopt_date,
    W,
    check.names = FALSE
  )
  rownames(W_correct) <- as.character(1:10)
  X_correct <- data.frame(
    adopt_date = adopt_date,
    popwt = 21:30,
    contr_cov_1 = c(rep(1:2, each = 3), rep(1, 4)),
    treat_cov_1 = c(rep(1, 6), 1, 2, 2, 1)
  )

  expect_s3_class(panel$Y, "data.frame")
  expect_s3_class(panel$W, "data.frame")
  expect_s3_class(panel$X, "data.frame")
  expect_equal(panel$Y, Y_correct)
  expect_equal(panel$W, W_correct)
  expect_equal(panel$X, X_correct)
})
