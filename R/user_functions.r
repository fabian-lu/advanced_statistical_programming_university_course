#' Fits a Local Maximum Likelihood (LML) regression model using B-splines
#'
#' This function fits a Local Maximum Likelihood (LML) regression model to the given data using B-splines.
#' It estimates the regression function by iteratively updating the mean and variance parameters, along with
#' the B-spline coefficients, until convergence is reached. The B-spline basis functions are used to model the
#' non-linear relationship between the predictors and the response.
#'
#' @param vector_x A vector of predictor values.
#' @param response_vector_y A vector of response values.
#' @param knot_count The number of knots to use in the B-spline basis.
#' @param poly_degree The degree of the B-spline polynomials.
#' @param order_difference_matrix_r The order of the difference matrix used for penalty.
#' @param max_iterations The maximum number of iterations for the algorithm.
#' @param lambda_mu The tuning parameter (lambda) for updating the mean parameter.
#' @param lambda_sigma The tuning parameter (lambda) for updating the variance parameter.
#' @param opt_method The optimization method to use for tuning parameter selection ('GCV' or 'GAIC').
#' @param lambda_acc The tuning parameter for the optimization method.
#'
#' @return A list containing the estimated mean (mu_hat), estimated variance (sigma_hat),
#' estimated B-spline coefficients (beta_hat), and estimated gamma coefficients for variance (gamma_hat).
#'
#' @examples
#' # Generate testing data
#' n <- 100
#' degree <- 3
#' noise_sd <- 0.4
#' number_knots <- 10
#' lower_boundary <- -5
#' upper_boundary <- 5
#'
#' data <- generate_polynomial_data(n, degree, noise_sd, lower_boundary, upper_boundary)
#'
#' # Fit the lmls_bspline model
#' result <- lmls_bspline(
#'   data$x, data$y, knot_count = 40, poly_degree = 3,
#'   order_difference_matrix_r = 2, max_iterations = 50, lambda_mu = 3, lambda_sigma = 500,
#'   opt_method = "GCV", lambda_acc = 50
#' )
#'
#' # Print the result
#' print(result)
#'
#' @export
#'
lmls_bspline <- function(
    vector_x, response_vector_y, knot_count = 40, poly_degree = 3,
    order_difference_matrix_r = 2, max_iterations = 50, lambda_mu = 3, lambda_sigma = 500,
    opt_method = "GCV", lambda_acc = 50
) {

  # Start with initialization
  init_result <- init(vector_x, response_vector_y, knot_count, poly_degree,
                      lambda_mu, lambda_sigma, order_difference_matrix_r)
  mu_hat <- init_result$mu_hat
  sigma_hat <- init_result$sigma_hat
  beta_hat <- init_result$beta_hat
  gamma_hat <- init_result$gamma_hat
  design_matrix_Z <- init_result$design_matrix_Z
  penalty_matrix_K <- init_result$penalty_matrix_K

  # Return rs algo outcome
  return(
    outer(
      mu_hat, sigma_hat, design_matrix_Z, vector_x, response_vector_y,
      beta_hat, gamma_hat, lambda_mu, lambda_sigma, max_iterations, penalty_matrix_K,
      penalty_matrix_K, knot_count, poly_degree, opt_method, lambda_acc
    )
  )
}

#' Plotting function for lmls_bspline results
#'
#' This function plots the estimated mean function along with confidence intervals
#' for the lmls_bspline results.
#'
#' @param result An object of class "lmls_bspline" containing the result of the lmls_bspline function.
#' @param data The original data used in the lmls_bspline function.
#'
#' @return None
#' @export
#'
#' @examples
#' data <- generate_polynomial_data(100, 3, 0.4, -5, 5)
#' result <- lmls_bspline(data$x, data$y)
#' plot_lmls_bspline(result, data)
#' @importFrom graphics lines points polygon
#' @importFrom grDevices rgb
plot_lmls_bspline <- function(result, data) {
  mu_hat <- result$mu_hat
  sigma_hat <- result$sigma_hat

  lower_bound <- mu_hat - 2 * sigma_hat
  upper_bound <- mu_hat + 2 * sigma_hat

  # Plotting the data points
  plot(data$x, data$y, pch = 16, main = "Estimated Mean with Confidence Intervals",
       xlab = "X", ylab = "Y")

  # Plotting the estimated mean function
  lines(data$x, mu_hat, type = "l", col = "blue", lwd = 3)

  # Plotting the confidence intervals
  lines(data$x, lower_bound, col = "green", lty = 2)  # Lower boundary
  lines(data$x, upper_bound, col = "red", lty = 2)    # Upper boundary

  # Adjusted polygon code with transparency
  polygon(c(data$x, rev(data$x)), c(lower_bound, rev(upper_bound)),
          col = rgb(0.5, 0.5, 0.5, alpha = 0.5), border = NA)

  # Overlaying the line and dots again to ensure they are not obscured
  lines(data$x, mu_hat, type = "l", col = "blue", lwd = 3)
  points(data$x, data$y, pch = 16)
}

