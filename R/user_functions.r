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
    vector_x, response_vector_y, kFnot_count = 40, poly_degree = 3,
    order_difference_matrix_r = 2, max_iterations = 50, lambda_mu = 3, lambda_sigma = 500,
    opt_method = "GCV", lambda_acc = 50
) {
  
  # Check input for vector_x
  if (class(vector_x) != "numeric") {
    stop("Input Argument for 'vector_x' must be numeric")
  }
  
  # Check input for response_vector_y
  if (class(response_vector_y) != "numeric") {
    stop("Input Argument for 'response_vector_y' must be numeric")
  }
  
  # Check input for knot_count
  if (class(knot_count) != "numeric") {
    stop("Input Argument for 'knot_count' must be numeric")
  }
  
  # Check input for poly_degree
  if (class(poly_degree) != "numeric") {
    stop("Input Argument for 'poly_degree' must be numeric")
  }
  
  # Check input for order_difference_matrix_r
  if (class(order_difference_matrix_r) != "numeric") {
    stop("Input Argument for 'order_difference_matrix_r' must be numeric")
  }
  
  # Check input for max_iterations
  if (class(max_iterations) != "numeric") {
    stop("Input Argument for 'max_iterations' must be numeric")
  }
  
  # Check input for lambda_mu
  if (class(lambda_mu) != "numeric") {
    stop("Input Argument for 'lambda_mu' must be numeric")
  }
  
  # Check input for lambda_sigma
  if (class(lambda_sigma) != "numeric") {
    stop("Input Argument for 'lambda_sigma' must be numeric")
  }
  
  # Check input for lambda_acc
  if (class(lambda_acc) != "numeric") {
    stop("Input Argument for 'lambda_acc' must be numeric")
  }

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

#' Plotting function for lmls_bspline results using ggplot2
#'
#' This function plots the estimated mean function along with confidence intervals
#' for the lmls_bspline results using ggplot2.
#'
#' @param result An object of class "lmls_bspline" containing the result of the lmls_bspline function.
#' @param data The original data used in the lmls_bspline function.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' data <- generate_polynomial_data(100, 3, 0.4, -5, 5)
#' result <- lmls_bspline(data$x, data$y)
#' plot_lmls_bspline(result, data)
#' @importFrom ggplot2 ggplot geom_point geom_line geom_ribbon labs theme_minimal
plot_lmls_bspline <- function(result, data) {
  # Input validation check
  if (class(result) != "pSplineLocationScale") {
    stop("The 'result' object needs to be of class pSplineLocationScale.")
  }
  
  # Extract relevant data
  mu_hat <- result$mu_hat
  sigma_hat <- result$sigma_hat
  x <- data$x
  y <- data$y
  
  # Calculate lower and upper bounds for confidence intervals
  alpha <- 0.05
  lower_bound <- mu_hat - qnorm(1 - alpha / 2) * sigma_hat
  upper_bound <- mu_hat + qnorm(1 - alpha / 2) * sigma_hat
  
  # Create a ggplot object
  gg <- ggplot() +
    geom_point(data = data, aes(x = x, y = y), size = 2, shape = 16) +
    geom_line(data = data, aes(x = x, y = mu_hat), color = "blue", size = 1) +
    geom_ribbon(data = data, aes(x = x, ymin = lower_bound, ymax = upper_bound),
                fill = "gray70", alpha = 0.5) +
    labs(
      title = "Estimated Mean with Confidence Intervals",
      x = "X",
      y = "Y"
    ) +
    theme_minimal()
  
  return(gg)
}

#' Summary function for lmls_bspline results
#'
#' This function provides a summary of the results obtained from the lmls_bspline function.
#'
#' @param lmls_object An object of class "lmls_bspline" obtained from lmls_bspline function.
#' @param digits Number of digits to round numerical values in the summary.
#'
#' @return A summary of the lmls_bspline object.
#' @export
#'
#' @examples
#' data <- generate_polynomial_data(100, 3, 0.4, -5, 5)
#' result <- lmls_bspline(data$x, data$y)
#' summary_lmls_bspline(result)
summary_lmls_bspline <- function(lmls_object, digits = 3) {
  if (class(lmls_object) != "list" || !all(names(lmls_object) %in% c("mu_hat", "sigma_hat", "beta_hat", "gamma_hat"))) {
    stop("lmls_object should be a list with components 'mu_hat', 'sigma_hat', 'beta_hat', and 'gamma_hat'")
  }
  
  cat("Local Maximum Likelihood (LML) Regression Summary\n")
  cat("------------------------------------------------\n")
  
  cat("\nOptimization Method: ", lmls_object$opt_method, "\n")
  
  cat("\nOptimal Lambda for Mu: ", format(signif(lmls_object$lambda_mu, digits)), "\n")
  
  cat("\nOptimal Lambda for Sigma: ", format(signif(lmls_object$lambda_sigma, digits)), "\n")
  
  cat("\nNumber of Knots: ", length(lmls_object$beta_hat) - 2, "\n")
  
  cat("\nDegree of B-spline Polynomials: ", lmls_object$poly_degree, "\n")
  
  cat("\nGlobal Deviance at Convergence: ", format(signif(lmls_object$GD_conv, digits)), "\n")
  
  cat("\nEstimated Parameters:\n")
  cat("Estimated Mean (mu_hat):\n")
  print(lmls_object$mu_hat, digits = digits)
  
  cat("\nEstimated Variance (sigma_hat):\n")
  print(lmls_object$sigma_hat, digits = digits)
  
  cat("\nEstimated B-spline Coefficients (beta_hat):\n")
  print(lmls_object$beta_hat, digits = digits)
  
  cat("\nEstimated Gamma Coefficients for Variance (gamma_hat):\n")
  print(lmls_object$gamma_hat, digits = digits)
}


