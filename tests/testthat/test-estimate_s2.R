## No Fixed and no Interactive Fixed Effects

test_that("OLS: No FE, no interactive FE. OK", {
  set.seed(42)
  N <- 100; T <- 100
  Y <- matrix(rnorm(N*T), nrow = N, ncol = T)
  W <- matrix(rbinom(n = N*T, size = 1, prob = 0.5), nrow = N, ncol = T)
  order_data <- order(rowSums(W), decreasing = TRUE)
  Y <- Y[order_data, ]; W <- W[order_data, ]

  s2 <- estimate_s2(Y, W)
  cond <- (s2 <= 1.1) && (s2 >= 0.9)
  expect_true(cond)
})

test_that("Weighted: No FE, no interactive FE. OK", {
  set.seed(142)
  N <- 100; T <- 100
  popwt <- rexp(N)
  Y <- matrix(rnorm(N*T), nrow = N, ncol = T)
  W <- matrix(rbinom(n = N*T, size = 1, prob = 0.5), nrow = N, ncol = T)
  popwt <- rexp(N)
  order_data <- order(rowSums(W), decreasing = TRUE)
  Y <- Y[order_data, ]; W <- W[order_data, ]; popwt <- popwt[order_data]

  s2 <- estimate_s2(Y, W, popwt)
  cond <- (s2 <= 1.1) && (s2 >= 0.9)
  expect_true(cond)
})

test_that("WLS: No FE, no interactive FE. OK", {
  set.seed(242)
  N <- 100; T <- 100
  Y <- matrix(rnorm(N*T), nrow = N, ncol = T)
  W <- matrix(rbinom(n = N*T, size = 1, prob = 0.5), nrow = N, ncol = T)
  popwt <- rexp(N)
  order_data <- order(rowSums(W), decreasing = TRUE)
  Y <- Y[order_data, ]; W <- W[order_data, ]; popwt <- popwt[order_data]

  s2 <- estimate_s2(Y, W, popwt, TRUE)
  cond <- (s2 <= 1.1) && (s2 >= 0.9)
  expect_true(cond)
})

## Fixed effects are present and Interactive Fixed Effects are no

test_that("OLS: FE, no interactive FE. OK", {
  set.seed(342)
  N <- 100; T <- 100
  FE <- outer(rexp(N, 1/7), rep(1, T)) + outer(rep(1, N), rexp(T, 1/5))
  Y <- matrix(rnorm(N * T), nrow = N, ncol = T) + FE
  W <- matrix(rbinom(n = N*T, size = 1, prob = 0.5), nrow = N, ncol = T)
  order_data <- order(rowSums(W), decreasing = TRUE)
  Y <- Y[order_data, ]; W <- W[order_data, ]

  s2 <- estimate_s2(Y, W)
  cond <- (s2 <= 1.1) && (s2 >= 0.9)
  expect_true(cond)
})

test_that("Weighted: No FE, no interactive FE is OK", {
  set.seed(142)
  N <- 100; T <- 100
  FE <- outer(rexp(N, 1/7), rep(1, T)) + outer(rep(1, N), rexp(T, 1/5))
  Y <- matrix(rnorm(N * T), nrow = N, ncol = T) + FE
  W <- matrix(rbinom(n = N*T, size = 1, prob = 0.5), nrow = N, ncol = T)
  popwt <- rexp(N)
  order_data <- order(rowSums(W), decreasing = TRUE)
  Y <- Y[order_data, ]; W <- W[order_data, ]; popwt <- popwt[order_data]

  s2 <- estimate_s2(Y, W, popwt)
  cond <- (s2 <= 1.1) && (s2 >= 0.9)
  expect_true(cond)
})

test_that("WLS: No FE, no interactive FE is OK", {
  set.seed(242)
  N <- 100; T <- 100
  FE <- outer(rexp(N, 1/7), rep(1, T)) + outer(rep(1, N), rexp(T, 1/5))
  Y <- matrix(rnorm(N * T), nrow = N, ncol = T) + FE
  W <- matrix(rbinom(n = N*T, size = 1, prob = 0.5), nrow = N, ncol = T)
  popwt <- rexp(N)
  order_data <- order(rowSums(W), decreasing = TRUE)
  Y <- Y[order_data, ]; W <- W[order_data, ]; popwt <- popwt[order_data]

  s2 <- estimate_s2(Y, W, popwt, TRUE)
  cond <- (s2 <= 1.1) && (s2 >= 0.9)
  expect_true(cond)
})




