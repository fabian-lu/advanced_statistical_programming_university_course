#' Creating Penalized B-Splines
#'
#' This function calculates the knots of the splines using a set of equidistant knots.
#' The number of inner knots is determined by knot_count, while the number of outer knots
#' is twice the degree of the polynomial.
#'
#' @param knot_count The number of knots (including inner and outer knots).
#' @param lower_bound The lower bound of the spline range.
#' @param upper_bound The upper bound of the spline range.
#' @param poly_degree The degree of the polynomial.
#'
#' @return The sequence of knots for the splines.
#'
create_knots_for_cox_de_boor <- function(
    knot_count, lower_bound, upper_bound, poly_degree
) {
  distance <- (upper_bound - lower_bound) / (knot_count - 1)
  inner_knots <- seq(from = lower_bound, to = upper_bound, length.out = knot_count)

  outer_knots_left <- seq(
    from = lower_bound - distance * poly_degree, by = distance,
    length.out = poly_degree
  )
  outer_knots_right <- seq(
    from = upper_bound + distance, by = distance, length.out = poly_degree
  )

  knots_sequence <- c(outer_knots_left, inner_knots, outer_knots_right)

  return(knots_sequence)
}

#' Indicator Function
#'
#' This function returns an indicator value based on the position relative to lower and upper bounds.
#' It returns 1 if the position is greater than or equal to the lower bound and less than the upper bound, and 0 otherwise.
#'
#' @param position The position value to evaluate.
#' @param lower The lower bound value.
#' @param upper The upper bound value.
#'
#' @return An indicator value (1 or 0) based on the position relative to the lower and upper bounds.
#'
indicator <- function(position, lower, upper) {
  return(as.integer(position >= lower & position < upper))
}



#' Cox de Boor Algorithm
#'
#' This function implements the Cox de Boor algorithm for calculating B-spline basis functions.
#' It computes the value of the B-spline basis function at a given position within a specific
#' range defined by lower and upper.
#'
#' @param position The position at which to evaluate the basis function.
#' @param knots The sequence of knots for the splines.
#' @param lower The lower index of the basis function.
#' @param upper The upper index of the basis function.
#'
#' @return The value of the B-spline basis function at the given position.
#'
cox_de_boor <- function(position, knots, lower, upper) {
  if (upper == 0) {
    return(indicator(position, knots[lower], knots[lower+1]))
  } else if (upper == 1) {
    term1 <- (position - knots[lower-1]) / (knots[lower] - knots[lower-1]) * indicator(position, knots[lower-1], knots[lower])
    term2 <- (knots[lower+1] - position) / (knots[lower+1] - knots[lower]) * indicator(position, knots[lower], knots[lower+1])
    return(term1 + term2)
  } else {
    term1 <- (position - knots[lower-upper]) / (knots[lower] - knots[lower-upper]) * cox_de_boor(position, knots, lower-1, upper-1)
    term2 <- (knots[lower+1] - position) / (knots[lower+1] - knots[lower+1-upper]) * cox_de_boor(position, knots, lower, upper-1)
    return(term1 + term2)
  }
}


#' B-Spline Basis Function Vector
#'
#' This function computes the B-spline basis function vector for a given set of knots and
#' a vector of x values.
#'
#' @param knots_sequence The sequence of knots for the splines.
#' @param vector_x The input vector of x values.
#' @param poly_degree The degree of the polynomial for the B-splines.
#' @param basis_function_knot_position The position of the basis function within the knot sequence.
#'
#' @return The B-spline basis function vector.
#'
create_basis_function_vector <- function(
    knots_sequence, vector_x, poly_degree, basis_function_knot_position
) {
  basis <- sapply(vector_x, function(x) cox_de_boor(x, knots_sequence, basis_function_knot_position, poly_degree))
  return(basis)
}


#' B-Spline Basis Function Matrix
#'
#' This function combines the B-spline basis function vectors into a matrix.
#' Each column of the matrix corresponds to a different basis function.
#' The resulting matrix is the design matrix Z for the B-splines.
#'
#' @param knot_count The number of knots for the splines.
#' @param lower_bound The lower bound of the spline range.
#' @param upper_bound The upper bound of the spline range.
#' @param poly_degree The degree of the polynomial for the B-splines.
#' @param vector_x The input vector of x values.
#'
#' @return The B-spline basis function matrix (design matrix Z).
#'
create_basis_function_matrix <- function(
    knot_count, lower_bound, upper_bound, poly_degree, vector_x
) {
  knots_sequence <- create_knots_for_cox_de_boor(knot_count, lower_bound, upper_bound, poly_degree)

  list_with_basis_functions <- list()
  start <- (1 + poly_degree)
  end <- start + (knot_count + poly_degree - 1) - 1
  for (i in start:end) {
    basis_function <- create_basis_function_vector(knots_sequence, vector_x, poly_degree, i)
    list_with_basis_functions[[i]] <- basis_function
  }

  vector_with_basis_functions <- do.call(cbind, list_with_basis_functions)
  vector_with_basis_functions[is.na(vector_with_basis_functions)] <- 0

  return(vector_with_basis_functions)
}

#' Difference Matrix
#'
#' This function creates the difference matrix used in penalized B-splines.
#' The order of the difference matrix determines the number of differences between
#' adjacent B-spline basis functions.
#'
#' @param order_of_difference_matrix The order of the difference matrix.
#' @param poly_degree The degree of the polynomial for the B-splines.
#' @param knot_count The number of knots for the splines.
#'
#' @return The difference matrix.
#'
create_diff_matrix <- function(
    order_of_difference_matrix, poly_degree, knot_count
) {
  if (order_of_difference_matrix == 1) {
    return(create_first_order_diff_matrix(poly_degree + knot_count - 1 - order_of_difference_matrix + 1))
  } else {
    return(create_first_order_diff_matrix(poly_degree + knot_count - 1 - order_of_difference_matrix + 1) %*%
             create_diff_matrix(order_of_difference_matrix - 1, poly_degree, knot_count))
  }
}

#' First Order Difference Matrix
#'
#' This function creates the first order difference matrix used in penalized B-splines.
#' It is a square matrix with dimensions (poly_degree + knot_count - 2) x (poly_degree + knot_count - 1).
#' The matrix contains -1 and 1 along the diagonal and zeros elsewhere.
#'
#' @param dimensions The dimensions of the difference matrix.
#'
#' @return The first order difference matrix.
#'
create_first_order_diff_matrix <- function(dimensions) {
  rows <- dimensions - 1
  cols <- dimensions

  difference_matrix <- matrix(0, nrow = rows, ncol = cols)

  offset <- 0
  for (row in 1:rows) {
    difference_matrix[row, offset + 1] <- -1
    difference_matrix[row, offset + 2] <- 1
    offset <- offset + 1
  }

  return(difference_matrix)
}

#' Penalty Matrix K
#'
#' This function computes the penalty matrix K used in penalized B-splines.
#' The penalty matrix is calculated by multiplying the transpose of the difference matrix
#' with the difference matrix itself.
#'
#' @param order_of_difference_matrix The order of the difference matrix.
#' @param poly_degree The degree of the polynomial for the B-splines.
#' @param knot_count The number of knots for the splines.
#'
#' @return The penalty matrix K.
#'
create_penalty_matrix_K <- function(
    order_of_difference_matrix, poly_degree, knot_count
) {
  difference_matrix <- create_diff_matrix(order_of_difference_matrix, poly_degree, knot_count)

  penalty_matrix_K <- t(difference_matrix) %*% difference_matrix

  return(penalty_matrix_K)
}

#' Penalized Least Squares B-Splines
#'
#' This function performs the penalized least squares regression using B-splines.
#' It computes the B-spline design matrix Z, the penalty matrix K, and the coefficients gamma_hat.
#' The hyperparameter lambda determines the amount of penalization.
#'
#' @param knot_count The number of knots for the splines.
#' @param lower_bound The lower bound of the spline range.
#' @param upper_bound The upper bound of the spline range.
#' @param poly_degree The degree of the polynomial for the B-splines.
#' @param vector_x The input vector of x values.
#' @param response_vector_y The response vector of y values.
#' @param hyperparameter_lambda The hyperparameter for penalization.
#'
#' @return A list containing the B-spline design matrix Z and the coefficients gamma_hat.
#'

PLS_b_Splines <- function(
    knot_count, lower_bound, upper_bound, poly_degree,
    vector_x, response_vector_y, hyperparameter_lambda
) {
  Z <- create_basis_function_matrix(
    knot_count, lower_bound, upper_bound, poly_degree,
    vector_x
  )
  K <- create_penalty_matrix_K(2, poly_degree, knot_count)
  gamma_hat <- solve(t(Z) %*% Z + hyperparameter_lambda * K) %*% t(Z) %*% response_vector_y

  result <- list(Z = Z, gamma_hat = gamma_hat)
  return(result)
}
