# Test for init function
test_that("init function initializes correctly", {
  vector_x <- c(1, 2, 3)
  response_vector_y <- c(1, 2, 3)
  result <- init(vector_x, response_vector_y, 3, 2, 0.5, 0.5, 2)
  expect_equal(length(result$mu_hat), 3)
  expect_equal(length(result$sigma_hat), 3)
  expect_equal(length(result$beta_hat), 4)
  expect_equal(length(result$gamma_hat), 4)
  expect_equal(dim(result$design_matrix_Z), c(3, 4))
})

# Test for global_deviance function
test_that("global_deviance calculates correctly", {
  response_vector_y <- c(1, 2, 3)
  mu_hat <- c(1, 2, 3)
  sigma_hat <- c(1, 1, 1)
  expect_equal(round(global_deviance(response_vector_y, mu_hat, sigma_hat),1), 5.5)
})

# Test for compute_lambda_gcv function
test_that("compute_lambda_gcv returns correct lambda", {
  design_matrix_Z <- matrix(c(1, 2, 3, 4.2, 5.3, 6, 7, 8, 9), ncol = 3)
  W <- diag(c(1, 1, 1))
  w <- c(1, 1, 1)
  z <- c(1, 2, 3)
  K <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), ncol = 3)
  lambda_acc <- 10
  n <- 3
  lambda_mu_flag <- TRUE
  epsilon <- 1e-5
  I <- diag(3)

  result <- compute_lambda_gcv(lambda_mu_flag, design_matrix_Z, W, w, z, K + epsilon * I, lambda_acc, n)
  expect_true(result >= 0 && result <= 50)
})


# Test for inner_sigma function
test_that("inner_sigma updates sigma_hat and gamma_hat", {
  mu_hat <- c(1, 2, 3)
  sigma_hat <- c(1, 1, 1)
  design_matrix_Z <- matrix(c(1, 2, 3, 4.2, 5.3, 6, 7, 8, 9), ncol = 3)
  response_vector_y <- c(1, 2, 3)
  gamma_hat <- c(1, 1)
  lambda_sigma <- 0.5
  K_sigma <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), ncol = 3)
  lambda_acc <- 10
  result <- inner_sigma(mu_hat, sigma_hat, design_matrix_Z, response_vector_y, gamma_hat, lambda_sigma, K_sigma, lambda_acc)
  expect_equal(length(result$sigma_hat), 3)
  expect_equal(length(result$gamma_hat), 3)
})


test_that("outer function updates mu_hat, sigma_hat, beta_hat, and gamma_hat", {
  mu_hat <- c(0.7, 1.8, 2.9)
  sigma_hat <- c(0.6, 0.7, 0.8)
  design_matrix_Z <- matrix(c(2, 3, 4, 5.2, 6.3, 7, 8, 9, 10), ncol = 3)
  vector_x <- c(1.1, 2.2, 3.3)
  response_vector_y <- c(1.2, 2.3, 3.4)
  beta_hat <- c(0.3, 0.4, 0.5)
  gamma_hat <- c(0.4, 0.5, 0.6)
  lambda_mu <- 0.6
  lambda_sigma <- 0.6
  max_iterations <- 7
  K_mu <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), ncol = 3)
  K_sigma <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), ncol = 3)
  knot_count <- 4
  poly_degree <- 3
  opt_method <- "BIC"
  lambda_acc <- 8
  result <- outer(mu_hat, sigma_hat, design_matrix_Z, vector_x, response_vector_y, beta_hat, gamma_hat, lambda_mu, lambda_sigma, max_iterations, K_mu, K_sigma, knot_count, poly_degree, opt_method, lambda_acc)
  expect_equal(length(result$mu_hat), 3)
  expect_equal(length(result$sigma_hat), 3)
  expect_equal(length(result$beta_hat), 3)
  expect_equal(length(result$gamma_hat), 3)
})

