test_that("cohort level from Y and X without covariates works fine", {

  adopt_date <- as.character(c(rep(2500, 5), rep(2014, 2), rep(2012, 3)))
  colnames <- as.character(2010:2019)
  Y <- matrix(1:100, nrow = 10, byrow = TRUE, dimnames = list(1:10, colnames))
  Y_df <- data.frame(adopt_date = adopt_date,  Y, check.names = FALSE)
  W <- matrix(
    c(
      rep(0, 50),
      rep(c(rep(0, 4), rep(1, 6)), 2),
      rep(c(rep(0, 2), rep(1, 8)), 3)
    ),
    nrow = 10, byrow = TRUE,
    dimnames = list(1:10, colnames)
  )

  W_df <- data.frame(adopt_date = adopt_date, W, check.names = FALSE)
  X <- data.frame(adopt_date = adopt_date, popwt = rep(1, 10))
  panel <- list(Y = Y_df, W = W_df, X = X, never_treat = "2500")
  panel_avg <- prepare_wide(panel, "cohort")

  Y_correct <- matrix(
    c(
      colSums(Y[1:5, ] / 5),
      colSums(Y[6:7, ] / 2),
      colSums(Y[8:10,] / 3)
    ),
    ncol = 10, byrow = TRUE,
    dimnames = list(c("2500", "2014", "2012"), colnames)
  )
  W_correct <- matrix(
    c(
      rep(0, 10),
      c(rep(0, 4), rep(1, 6)),
      c(rep(0, 2), rep(1, 8))
    ),
    ncol = 10, byrow = TRUE,
    dimnames = list(c("2500", "2014", "2012"), colnames)
  )
  coh_correct <- c(5, 2, 3)

  expect_true(is.matrix(panel_avg$Y_avg))
  expect_true(is.matrix(panel_avg$W_avg))
  expect_equal(panel_avg$Y_avg, Y_correct)
  expect_equal(panel_avg$W_avg, W_correct)
  expect_equal(panel_avg$coh_weights, coh_correct)
})

test_that("cohort level from Y_wt without covariates works fine", {

  adopt_date <- as.character(c(rep(2500, 5), rep(2014, 2), rep(2012, 3)))
  colnames <- as.character(2010:2019)
  Y <- matrix(1:100, nrow = 10, byrow = TRUE, dimnames = list(1:10, colnames))
  # This is the difference with the previous test
  Y_wt <- data.frame(
    adopt_date = adopt_date,
    popwt = rep(1, 10),
    Y,
    check.names = FALSE
  )
  W <- matrix(
    c(
      rep(0, 50),
      rep(c(rep(0, 4), rep(1, 6)), 2),
      rep(c(rep(0, 2), rep(1, 8)), 3)
    ),
    nrow = 10, byrow = TRUE,
    dimnames = list(1:10, colnames)
  )
  W_df <- data.frame(adopt_date = adopt_date, W, check.names = FALSE)

  panel <- list(Y_wt = Y_wt, W = W_df, never_treat = "2500")
  panel_avg <- prepare_wide(panel)

  Y_correct <- matrix(
    c(
      colSums(Y[1:5, ] / 5),
      colSums(Y[6:7, ] / 2),
      colSums(Y[8:10, ] / 3)
    ),
    ncol = 10, byrow = TRUE,
    dimnames = list(c("2500", "2014", "2012"), colnames)
  )
  W_correct <- matrix(
    c(
      rep(0, 10),
      c(rep(0, 4), rep(1, 6)),
      c(rep(0, 2), rep(1, 8))
    ),
    ncol = 10, byrow = TRUE,
    dimnames = list(c("2500", "2014", "2012"), colnames)
  )
  coh_correct <- c(5, 2, 3)

  expect_true(is.matrix(panel_avg$Y_avg))
  expect_true(is.matrix(panel_avg$W_avg))
  expect_equal(panel_avg$Y_avg, Y_correct)
  expect_equal(panel_avg$W_avg, W_correct)
  expect_equal(panel_avg$coh_weights, coh_correct)
})

test_that("cohort level from Y and X with control covariate works fine", {

  adopt_date <- as.character(c(rep(2500, 5), rep(2014, 2), rep(2012, 3)))
  colnames <- as.character(2010:2019)
  Y <- matrix(1:100, nrow = 10, byrow = TRUE, dimnames = list(1:10, colnames))
  Y_df <- data.frame(adopt_date = adopt_date,  Y, check.names = FALSE)
  W <- matrix(
    c(
      rep(0, 50),
      rep(c(rep(0, 4), rep(1, 6)), 2),
      rep(c(rep(0, 2), rep(1, 8)), 3)
    ),
    nrow = 10, byrow = TRUE,
    dimnames = list(1:10, colnames)
  )

  W_df <- data.frame(adopt_date = adopt_date, W, check.names = FALSE)
  X <- data.frame(
    adopt_date = adopt_date,
    contr_cov_1 = c(rep(1, 2), rep(2, 3), rep(1, 5)),
    popwt = rep(1, 10))
  panel <- list(Y = Y_df, W = W_df, X = X, never_treat = "2500")
  panel_avg <- prepare_wide(panel)
  Y_correct_1 <- matrix(
    c(
      colSums(Y[1:2, ] / 2),
      colSums(Y[3:5, ] / 3),
      colSums(Y[6:7, ] / 2),
      colSums(Y[8:10,] / 3)
    ),
    ncol = 10, byrow = TRUE,
    dimnames = list(c("2500_1", "2500_2", "2014", "2012"), colnames)
  )
  Y_correct_2 <- matrix(
    c(
      colSums(Y[3:5, ] / 3),
      colSums(Y[1:2, ] / 2),
      colSums(Y[6:7, ] / 2),
      colSums(Y[8:10,] / 3)
    ),
    ncol = 10, byrow = TRUE,
    dimnames = list(c("2500_2", "2500_1", "2014", "2012"), colnames)
  )
  W_correct_1 <- matrix(
    c(
      rep(0, 10),
      rep(0, 10),
      c(rep(0, 4), rep(1, 6)),
      c(rep(0, 2), rep(1, 8))
    ),
    ncol = 10, byrow = TRUE,
    dimnames = list(c("2500_1", "2500_2", "2014", "2012"), colnames)
  )
  W_correct_2 <- matrix(
    c(
      rep(0, 10),
      rep(0, 10),
      c(rep(0, 4), rep(1, 6)),
      c(rep(0, 2), rep(1, 8))
    ),
    ncol = 10, byrow = TRUE,
    dimnames = list(c("2500_2", "2500_1", "2014", "2012"), colnames)
  )
  coh_correct_1 <- c(2, 3, 2, 3)
  coh_correct_2 <- c(3, 2, 2, 3)

  expect_true(is.matrix(panel_avg$Y_avg))
  expect_true(is.matrix(panel_avg$W_avg))

  #Y_avg
  expect_true(
    any(c(identical(panel_avg$Y_avg, Y_correct_1),
          identical(panel_avg$Y_avg, Y_correct_2)))
  )
  #W_avg
  expect_true(
    any(c(identical(panel_avg$W_avg, W_correct_1),
          identical(panel_avg$W_avg, W_correct_2)))
  )
  expect_true(
    any(c(identical(panel_avg$coh_weights, coh_correct_1),
          identical(panel_avg$coh_weights, coh_correct_2)))
  )
})

test_that("cohort level from Y_wt with control covariate works fine", {

  adopt_date <- as.character(c(rep(2500, 5), rep(2014, 2), rep(2012, 3)))
  colnames <- as.character(2010:2019)
  Y <- matrix(1:100, nrow = 10, byrow = TRUE, dimnames = list(1:10, colnames))
  Y_wt <- data.frame(
    adopt_date = adopt_date,
    popwt = rep(1, 10),
    contr_cov_1 = c(rep(1, 2), rep(2, 3), rep(1, 5)),
    Y,
    check.names = FALSE
  )
  W <- matrix(
    c(
      rep(0, 50),
      rep(c(rep(0, 4), rep(1, 6)), 2),
      rep(c(rep(0, 2), rep(1, 8)), 3)
    ),
    nrow = 10, byrow = TRUE,
    dimnames = list(1:10, colnames)
  )

  W_df <- data.frame(adopt_date = adopt_date, W, check.names = FALSE)
  panel <- list(Y_wt = Y_wt, W = W_df, never_treat = "2500")
  panel_avg <- prepare_wide(panel)

  Y_correct_1 <- matrix(
    c(
      colSums(Y[1:2, ] / 2),
      colSums(Y[3:5, ] / 3),
      colSums(Y[6:7, ] / 2),
      colSums(Y[8:10,] / 3)
    ),
    ncol = 10, byrow = TRUE,
    dimnames = list(c("2500_1", "2500_2", "2014", "2012"), colnames)
  )
  Y_correct_2 <- matrix(
    c(
      colSums(Y[3:5, ] / 3),
      colSums(Y[1:2, ] / 2),
      colSums(Y[6:7, ] / 2),
      colSums(Y[8:10,] / 3)
    ),
    ncol = 10, byrow = TRUE,
    dimnames = list(c("2500_2", "2500_1", "2014", "2012"), colnames)
  )
  W_correct_1 <- matrix(
    c(
      rep(0, 10),
      rep(0, 10),
      c(rep(0, 4), rep(1, 6)),
      c(rep(0, 2), rep(1, 8))
    ),
    ncol = 10, byrow = TRUE,
    dimnames = list(c("2500_1", "2500_2", "2014", "2012"), colnames)
  )
  W_correct_2 <- matrix(
    c(
      rep(0, 10),
      rep(0, 10),
      c(rep(0, 4), rep(1, 6)),
      c(rep(0, 2), rep(1, 8))
    ),
    ncol = 10, byrow = TRUE,
    dimnames = list(c("2500_2", "2500_1", "2014", "2012"), colnames)
  )
  coh_correct_1 <- c(2, 3, 2, 3)
  coh_correct_2 <- c(3, 2, 2, 3)

  expect_true(is.matrix(panel_avg$Y_avg))
  expect_true(is.matrix(panel_avg$W_avg))
  #Y_avg
  expect_true(
    any(c(identical(panel_avg$Y_avg, Y_correct_1),
          identical(panel_avg$Y_avg, Y_correct_2)))
  )
  #W_avg
  expect_true(
    any(c(identical(panel_avg$W_avg, W_correct_1),
          identical(panel_avg$W_avg, W_correct_2)))
  )
  #coh
  expect_true(
    any(c(identical(panel_avg$coh_weights, coh_correct_1),
          identical(panel_avg$coh_weights, coh_correct_2)))
  )
})

test_that("cohort level from Y with X with contr and treat covs and non unit weights (population) works fine", {
  adopt_date <- as.character(c(rep(2500, 5), rep(2014, 5), rep(2012, 10)))
  colnames <- as.character(2010:2019)
  Y <- matrix(1:200, nrow = 20, byrow = TRUE, dimnames = list(1:20, colnames))
  Y_df <- data.frame(adopt_date = adopt_date,  Y, check.names = FALSE)
  W <- matrix(
    c(
      rep(0, 50),
      rep(c(rep(0, 4), rep(1, 6)), 5),
      rep(c(rep(0, 2), rep(1, 8)), 10)
    ),
    nrow = 20, byrow = TRUE,
    dimnames = list(1:20, colnames)
  )

  W_df <- data.frame(adopt_date = adopt_date, W, check.names = FALSE)
  X <- data.frame(
    adopt_date = adopt_date,
    popwt = c(c(3, 1), rep(2, 3), rep(2, 3), rep(1, 2), rep(1, 10)),
    contr_cov_1 = c(rep(1, 2), rep(2, 3), rep(1, 15)),
    treat_cov_1 = c(rep(1, 5), rep(1, 5), rep(1, 4), rep(2, 6))
  )
  panel <- list(Y = Y_df, W = W_df, X = X, never_treat = "2500")
  panel_avg <- prepare_wide(panel, "cohort")

  Y_correct <- matrix(c(
      c(1:10 + 2.5),
      c(31:40),
      c(1:10 + 66.25),
      c(116:125),
      c(166:175)
    ),
    byrow = TRUE, ncol = 10,
    dimnames = list(c("2500_1", "2500_2", "2014_1", "2012_1", "2012_2"), colnames)
  )
  W_correct <- matrix(c(
      rep(0, 20),
      c(rep(0, 4), rep(1, 6)),
      rep(c(rep(0, 2), rep(1, 8)), 2)
    ),
    byrow = TRUE, ncol = 10,
    dimnames = list(c("2500_1", "2500_2", "2014_1", "2012_1", "2012_2"), colnames)
  )

  expect_true(is.matrix(panel_avg$Y_avg))
  expect_true(is.matrix(panel_avg$W_avg))
  expect_equal(panel_avg$Y_avg, Y_correct)
  expect_equal(panel_avg$W_avg, W_correct)
  expect_equal(panel_avg$coh_weights, c(4, 6, 8, 4, 6))
})

test_that("cohort level from Y_wt with contr and treat covs and non unit weights (population) works fine", {
  adopt_date <- as.character(c(rep(2500, 5), rep(2014, 5), rep(2012, 10)))
  colnames <- as.character(2010:2019)
  Y <- matrix(1:200, nrow = 20, byrow = TRUE, dimnames = list(1:20, colnames))
  Y_wt <- data.frame(
    adopt_date = adopt_date,
    Y,
    popwt = c(c(3, 1), rep(2, 3), rep(2, 3), rep(1, 2), rep(1, 10)),
    contr_cov_1 = c(rep(1, 2), rep(2, 3), rep(1, 15)),
    treat_cov_1 = c(rep(1, 5), rep(1, 5), rep(1, 4), rep(2, 6)),
    check.names = FALSE)
  W <- matrix(
    c(
      rep(0, 50),
      rep(c(rep(0, 4), rep(1, 6)), 5),
      rep(c(rep(0, 2), rep(1, 8)), 10)
    ),
    nrow = 20, byrow = TRUE,
    dimnames = list(1:20, colnames)
  )

  W_df <- data.frame(adopt_date = adopt_date, W, check.names = FALSE)
  panel <- list(Y_wt = Y_wt, W = W_df, never_treat = "2500")
  panel_avg <- prepare_wide(panel, "cohort")

  Y_correct <- matrix(c(
      c(1:10 + 2.5),
      c(31:40),
      c(1:10 + 66.25),
      c(116:125),
      c(166:175)
    ),
    byrow = TRUE, ncol = 10,
    dimnames = list(c("2500_1", "2500_2", "2014_1", "2012_1", "2012_2"), colnames)
  )
  W_correct <- matrix(c(
      rep(0, 20),
      c(rep(0, 4), rep(1, 6)),
      rep(c(rep(0, 2), rep(1, 8)), 2)
    ),
    byrow = TRUE, ncol = 10,
    dimnames = list(c("2500_1", "2500_2", "2014_1", "2012_1", "2012_2"), colnames)
  )

  expect_true(is.matrix(panel_avg$Y_avg))
  expect_true(is.matrix(panel_avg$W_avg))
  expect_equal(panel_avg$Y_avg, Y_correct)
  expect_equal(panel_avg$W_avg, W_correct)
  expect_equal(panel_avg$coh_weights, c(4, 6, 8, 4, 6))
})


# unit level
test_that("unit level from Y with X works fine", {
  adopt_date <- as.character(c(rep(2500, 5), rep(2014, 2), rep(2012, 3)))
  colnames <- as.character(2010:2019)
  Y <- matrix(1:100, nrow = 10, byrow = TRUE, dimnames = list(1:10, colnames))
  Y_df <- data.frame(adopt_date = adopt_date,  Y, check.names = FALSE)
  W <- matrix(
    c(
      rep(0, 50),
      rep(c(rep(0, 4), rep(1, 6)), 2),
      rep(c(rep(0, 2), rep(1, 8)), 3)
    ),
    nrow = 10, byrow = TRUE,
    dimnames = list(1:10, colnames)
  )

  W_df <- data.frame(adopt_date = adopt_date, W, check.names = FALSE)
  X <- data.frame(adopt_date = adopt_date, popwt = rep(1, 10))
  panel <- list(Y = Y_df, W = W_df, X = X, never_treat = "2500")
  panel_avg <- prepare_wide(panel, "unit")

  expect_true(is.matrix(panel_avg$Y_avg))
  expect_true(is.matrix(panel_avg$W_avg))
  expect_equal(panel_avg$Y_avg, Y)
  expect_equal(panel_avg$W_avg, W)
  expect_equal(panel_avg$coh_weights, rep(1, 10))
})

test_that("unit level from Y_wt works fine", {
  adopt_date <- as.character(c(rep(2500, 5), rep(2014, 2), rep(2012, 3)))
  colnames <- as.character(2010:2019)
  Y <- matrix(1:100, nrow = 10, byrow = TRUE, dimnames = list(1:10, colnames))
  Y_wt <- data.frame(
    adopt_date = adopt_date,
    Y,
    popwt = rep(1, 10),
    check.names = FALSE
  )
  W <- matrix(
    c(
      rep(0, 50),
      rep(c(rep(0, 4), rep(1, 6)), 2),
      rep(c(rep(0, 2), rep(1, 8)), 3)
    ),
    nrow = 10, byrow = TRUE,
    dimnames = list(1:10, colnames)
  )

  W_df <- data.frame(adopt_date = adopt_date, W, check.names = FALSE)
  panel <- list(Y_wt = Y_wt, W = W_df, never_treat = "2500")
  panel_avg <- prepare_wide(panel, "unit", FALSE)
  expect_true(is.matrix(panel_avg$Y_avg))
  expect_true(is.matrix(panel_avg$W_avg))
  expect_equal(panel_avg$Y_avg, Y)
  expect_equal(panel_avg$W_avg, W)
  expect_equal(panel_avg$coh_weights, rep(1, 10))
})
