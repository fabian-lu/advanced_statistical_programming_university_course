# Test for create_knots_for_cox_de_boor function
test_that("create_knots_for_cox_de_boor returns correct knots", {
  result <- create_knots_for_cox_de_boor(3, 0, 2, 2)
  expected <- c(-2, -1, 0, 1, 2, 3, 4)
  expect_equal(result, expected)
})

# Test for indicator function
test_that("indicator returns correct values", {
  expect_equal(indicator(5, 1, 10), 1)
  expect_equal(indicator(0, 1, 10), 0)
  expect_equal(indicator(10, 1, 10), 0)
})


# Test for create_basis_function_vector function
test_that("create_basis_function_vector returns correct values", {
  knots <- c(0, 1, 2, 3, 4)
  vector_x <- c(1, 1.5, 2)
  result <- create_basis_function_vector(knots, vector_x, 2, 3)
  expected <- c(0.5, 0.75, 0.5)
  expect_equal(result, expected)
})


# Test for create_first_order_diff_matrix function
test_that("create_first_order_diff_matrix returns correct matrix", {
  result <- create_first_order_diff_matrix(3)
  expected <- matrix(c(-1, 0, 1, -1, 0, 1), ncol = 3)
  expect_equal(result, expected)
})


# Test for create_penalty_matrix_K function
test_that("create_penalty_matrix_K returns correct matrix", {
  result <- create_penalty_matrix_K(2, 2, 3)
  expected <- matrix(c(1, -2, 1, 0, -2, 5, -4, 1, 1, -4, 5, -2, 0, 1, -2, 1), ncol = 4)
  expect_equal(result, expected)
})



# Test for PLS_b_Splines function
test_that("PLS_b_Splines returns correct values", {
  vector_x <- c(1, 1.5, 2)
  response_vector_y <- c(1, 2, 3)
  result <- PLS_b_Splines(3, 0, 2, 2, vector_x, response_vector_y, 0.5)
  expect_equal(dim(result$Z), c(3, 4))
  expect_equal(length(result$gamma_hat), 4)
})

