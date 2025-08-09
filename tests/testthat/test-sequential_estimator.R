# Check that there are necesssary outputs
test_that("All the essential elements are present in the output: 'sdid'.", {
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

  panel <- to_wide(
    panel_long,
    unit = "unit", time = "time",
    outcome = 1, treatment = 4,
    population = 5,
    contr_covs = "contr_cov_1",
    treat_covs = c("treat_cov_1"),
    sort = TRUE)

  panel_avg <- prepare_wide(panel, "cohort")
  est <- sequential_estimator(panel_avg, "cohort", 0.1)

  expect_s3_class(est, "sequential_estimate")
  expect_equal(attr(est, "type"), "sdid")
  expect_equal(attr(est, "level"), "cohort")
  expect_equal(attr(est, "s2"), 0.1)
})

test_that("All the essential elements are present in the output: 'did'.", {
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

  panel <- to_wide(
    panel_long,
    unit = "unit", time = "time",
    outcome = 1, treatment = 4,
    population = 5,
    contr_covs = "contr_cov_1",
    treat_covs = c("treat_cov_1"),
    sort = TRUE)

  panel_avg <- prepare_wide(panel, "cohort")
  est <- sequential_estimator(panel_avg, "cohort", type = "did")

  expect_s3_class(est, "sequential_estimate")
  expect_equal(attr(est, "type"), "did")
  expect_equal(attr(est, "level"), "cohort")
  expect_true(attr(est, "s2") < 1e-8)
})

test_that("All the essential elements are present in the output: 'both'.", {
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

  panel <- to_wide(
    panel_long,
    unit = "unit", time = "time",
    outcome = 1, treatment = 4,
    population = 5,
    contr_covs = "contr_cov_1",
    treat_covs = c("treat_cov_1"),
    sort = TRUE)

  panel_avg <- prepare_wide(panel, "unit")
  est <- sequential_estimator(panel_avg, "cohort", type = "did")

  expect_s3_class(est, "sequential_estimate")
  expect_equal(attr(est, "type"), "did")
  expect_equal(attr(est, "level"), "cohort")
  expect_true(attr(est, "s2") < 1e-8)
})

test_that("All the essential elements are present in the output on `unit` level and 'both' estimates.", {
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

  panel <- to_wide(
    panel_long,
    unit = "unit", time = "time",
    outcome = 1, treatment = 4,
    population = 5,
    contr_covs = "contr_cov_1",
    treat_covs = c("treat_cov_1"),
    sort = TRUE)

  panel_avg <- prepare_wide(panel, "unit")
  est <- sequential_estimator(panel_avg, "cohort", type = "both")

  expect_s3_class(est, "sequential_estimate")
  expect_equal(names(est), c("tau_sdid", "tau_did"))
  expect_equal(attr(est, "type"), "both")
  expect_equal(attr(est, "level"), "cohort")
  expect_true(attr(est, "s2") < 1e-8)

})

test_that("The estimated effect is around zero for FE and noise", {
  panel_long <- data.frame(
    outcome = 1:5000,
    unit = rep(1:500, each = 10),
    time = rep(1:10, 500),
    w = c(rep(0, 3000),
          rep(c(rep(0, 6), rep(1, 4)), 50),
          rep(c(rep(0, 4), rep(1, 6)), 100),
          rep(c(rep(0, 2), rep(1, 8)), 50)),
    population = rep(1:500, each = 10),
    contr_cov_1 = c(rep(1, 1500), rep(2, 1500), rep(1, 2000)),
    treat_cov_1 = c(rep(1, 3000), rep(1, 500), rep(2, 1000), rep(1, 500))
  )

  panel <- to_wide(
    panel_long,
    unit = "unit", time = "time",
    outcome = 1, treatment = 4,
    population = 5,
    contr_covs = "contr_cov_1",
    treat_covs = c("treat_cov_1"),
    sort = TRUE)

  panel_avg <- prepare_wide(panel, "unit")
  panel_avg$Y_avg <- panel_avg$Y_avg + matrix(rnorm(500*10), ncol = 10)

  est <- sequential_estimator(panel_avg, "cohort", type = "both")

  expect_s3_class(est, "sequential_estimate")
  expect_true(sum(est$tau_sdid^2)^0.5 < 1)
  expect_true(sum(est$tau_did^2)^0.5 < 1)
  expect_equal(names(est), c("tau_sdid", "tau_did"))
  expect_equal(attr(est, "type"), "both")
  expect_equal(attr(est, "level"), "cohort")
  expect_true((attr(est, "s2") < 1.1) & (attr(est, "s2") > 0.9))
})
