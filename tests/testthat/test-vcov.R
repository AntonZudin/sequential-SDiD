test_that("All essential elements are present", {
  set.seed(42)
  panel_long <- data.frame(
    outcome = 1:1000,
    unit = rep(1:100, each = 10),
    time = rep(1:10, 100),
    w = c(rep(0, 600),
          rep(c(rep(0, 6), rep(1, 4)), 10),
          rep(c(rep(0, 4), rep(1, 6)), 20),
          rep(c(rep(0, 2), rep(1, 8)), 10)),
    population = rep(1:100, each = 10),
    contr_cov_1 = c(rep(1, 300), rep(2, 300), rep(1, 400)),
    treat_cov_1 = c(rep(1, 600), rep(1, 100), rep(2, 200), rep(1, 100))
  )

  panel <- to_wide(
    panel_long,
    unit = "unit", time = "time",
    outcome = 1, treatment = 4,
    population = 5,
    contr_covs = "contr_cov_1",
    treat_covs = c("treat_cov_1"),
    sort = TRUE)

  panel$Y[, 2:11] <- panel$Y[, 2:11] + matrix(rnorm(100*10), ncol = 10)
  panel_avg <- prepare_wide(panel, "cohort")
  s2 <- estimate_s2(panel$Y, panel$W)
  est <- sequential_estimator(panel_avg, s2 = s2, level =  "cohort", type = "both")
  se <- vcov(est, panel, B = 1000)

  expect_equal(names(se), c("se_sdid", "se_did"))
  expect_equal(length(se$se_sdid), 8)
  expect_equal(length(se$se_did), 8)

  se <- vcov(est, panel, B = 1000, return_tau_b = TRUE)
  expect_equal(names(se), c("se", "tau_b"))
  expect_equal(names(se$se), c("se_sdid", "se_did"))
  expect_equal(names(se$tau_b), c("tau_b_sdid", "tau_b_did"))
  expect_equal(length(se$se$se_sdid), 8)
  expect_equal(length(se$se$se_did), 8)
  expect_equal(dim(se$tau_b$tau_b_sdid), c(8, 1000))
  expect_equal(dim(se$tau_b$tau_b_did), c(8, 1000))
})
