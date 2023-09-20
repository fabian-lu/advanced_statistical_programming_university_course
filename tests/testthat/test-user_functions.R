test_data <- generate_polynomial_data(100, 3, 0.4, -5, 5)

# Test for lmls_bspline
test_that("lmls_bspline returns correct object type", {
  result <- lmls_bspline(test_data$x, test_data$y)
  expect_s3_class(result, "pSplineLocationScale")
})

test_that("lmls_bspline returns object with expected components", {
  result <- lmls_bspline(test_data$x, test_data$y)
  expect_named(result, c("mu_hat", "sigma_hat", "beta_hat", "gamma_hat", "GD_conv", "iterations", "poly_degree", "max_iterations", "lambda_mu", "lambda_sigma", "opt_method"))
})


test_that("lmls_bspline throws error for incorrect input types", {
  expect_error(lmls_bspline(character(100), test_data$y), "Input Argument for 'vector_x' must be numeric")
})

# Test for plot_lmls_bspline
test_that("plot_lmls_bspline returns ggplot object", {
  result <- lmls_bspline(test_data$x, test_data$y)
  plot_result <- plot_lmls_bspline(result, test_data)
  expect_s3_class(plot_result, "ggplot")
})

test_that("plot_lmls_bspline throws error for incorrect input types", {
  expect_error(plot_lmls_bspline(data.frame(), test_data), "The 'result' object needs to be of class pSplineLocationScale.")
})

# Test for summary_lmls_bspline
test_that("summary_lmls_bspline prints expected components", {
  result <- lmls_bspline(test_data$x, test_data$y)
  captured_output <- capture.output(summary_lmls_bspline(result))
  # Join the captured output into a single string
  output_string <- paste(captured_output, collapse = "\n")
  # Check if the output string contains "Estimated Mean (mu_hat):"
  expect_true(grepl("Estimated Mean \\(mu_hat\\):", output_string))
})

test_that("summary_lmls_bspline throws error for incorrect input types", {
  expect_error(summary_lmls_bspline(data.frame()), "lmls_object should be a list with the required components")
})
